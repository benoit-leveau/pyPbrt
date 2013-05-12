"""Grid Accelerator."""

from core.pbrt import clamp, round_to_int
from core.rwlock import RWLock
from core.primitive import Aggregate
from core.geometry import BBox, Vector, union


class Voxel(object):

    """Class describing a Voxel."""

    def __init__(self, primitive=None):
        """Default constructor for Voxel."""
        self.primitives = []
        if primitive:
            self.primitives.append(primitive)
        self.all_can_intersect = False

    def size(self):
        """Return the number of primitives in the voxel."""
        return len(self.primitives)

    def add_primitive(self, primitive):
        """Add the primitive to the voxel."""
        self.primitives.append(primitive)
        
    def intersect(self, ray, intersection, lock):
        """Compute the intersection with the enclosed primitives."""
        if not self.all_can_intersect:
            self._fully_refine(lock)

        # loop over primitives in voxel and find closest (if any) intersection
        hit_something = False
        for primitive in self.primitives:
            if primitive.intersect(ray, intersection):
                # don't stop at first hit, we may have a nearer intersection
                hit_something = True

        # intersection should have been filled at this stage
        return hit_something

    def intersect_p(self, ray, lock):
        """Return True if the ray intersects any primitive in the voxel."""
        if not self.all_can_intersect:
            self._fully_refine(lock)
        
        for primitive in self.primitives:
            if primitive.intersect_p(ray):
                # stops at first intersection found
                return True

        # no intersection
        return False

    def _fully_refine(self, lock):
        """Fully refine all primitives in the voxel."""
        lock.promote()
        for i in range(len(self.primitives)):
            primitive = self.primitives[i]
            # refine primitive if it's not intersectable
            if not primitive.can_intersect():
                new_primitives = []
                primitive.fully_refine(new_primitives)
                if len(new_primitives) == 1:
                    # use the primitive as is in the voxel
                    self.primitives[i] = new_primitives[0]
                else:
                    # put the primitives in a new grid accel
                    self.primitives[i] = GridAccel(new_primitives, False)
        self.all_can_intersect = True
        lock.demote()
        
class GridAccel(Aggregate):

    """Class describing a GridAccel."""

    def __init__(self, primitives, refine_immediately):
        """Default constructor for GridAccel."""
        # initialize self.primitives with primitives for grid
        if refine_immediately:
            self.primitives = []
            for primitive in primitives:
                primitive.fully_refine(self.primitives)
        else:
            self.primitives = list(primitives)

        # compute bounds and choose grid resolution
        self.bounds = BBox()
        for primitive in self.primitives:
            self.bounds = union(self.bounds, primitive.world_bound())
        delta = self.bounds.p_max - self.bounds.p_min
        
        # find voxels_per_unit_dist for grid
        max_axis = self.bounds.maximum_extent()
        inv_max_width = 1.0 / delta[max_axis]
        cube_root = 3.0 * pow(len(self.primitives), 1.0/3.0)
        voxels_per_unit_dist = cube_root * inv_max_width
        self.n_voxels = []
        for axis in range(3):
            self.n_voxels.append(clamp(
                round_to_int(delta[axis] * voxels_per_unit_dist), 1, 64))

        # compute voxel widths and allocate voxels
        self.width = Vector()
        self.inv_width = Vector()
        for axis in range(3):
            self.width[axis] = delta[axis] / self.n_voxels[axis]
            if self.width[axis] == 0.0:
                self.inv_width[axis] = 0.0
            else:
                self.inv_width[axis] = 1.0 / self.width[axis]
        nv = self.n_voxels[0] * self.n_voxels[1] * self.n_voxels[2]

        # array of voxels, initialized at None
        self.voxels = [None] * nv
        
        # add primitives to grid voxels
        for primitive in self.primitives:
            # find voxel extent of primitive
            primitive_bound = primitive.world_bound()
            v_min = []
            v_max = []
            for axis in range(3):
                v_min.append(self._pos_to_voxel(primitive_bound.p_min, axis))
                v_max.append(self._pos_to_voxel(primitive_bound.p_max, axis))

            # add primitive to overlapping voxels
            for z in range(v_min[2], v_max[2]):
                for y in range(v_min[1], v_max[1]):
                    for x in range(v_min[0], v_max[0]):
                        index = self._offset(x, y, z)
                        if self.voxels[index] is None:
                            self.voxels[index] = Voxel(primitive)
                        else:
                            self.voxels[index].add_primitive(primitive)

        # create reader-writer mutex for grid
        self.rw_lock = RWLock()
        
    def world_bound(self):
        """Return the bounding box in world space."""
        return self.bounds

    def can_intersect(self):
        """Return True if the aggregate can intersect."""
        return True

    def intersect(self, ray, intersection):
        """Compute the intersection with the primitives."""
        # check ray against overall grid bounds
        if self.bounds.inside(ray(ray.mint)):
            ray_t = ray.mint
        else:
            intersected, t0, t1 = self.bounds.intersect_p(ray)
            if not intersected:
                self.rw_lock.release()
                return False
            ray_t = t0

        grid_intersect = ray(ray_t)

        # set up 3D DDA (Digital Differential Analyzer) for ray
        pos = []
        next_crossing_t = []
        delta_t = []
        step = []
        out = []
        for axis in range(3):
            # compute current voxel for axis
            pos.append(self._pos_to_voxel(grid_intersect, axis))
            if ray.d[axis] == 0.0:
                next_crossing_t.append(float('inf'))
                delta_t.append(float('inf'))
                step.append(1)
                out.append(self.n_voxels[axis])
            elif ray.d[axis] > 0.0:
                # handle ray with positive direction for voxel stepping
                next_crossing_t.append(
                    ray_t + (self._voxel_to_pos(pos[axis]+1, axis) - \
                             grid_intersect[axis]) / ray.d[axis])
                delta_t.append(self.width[axis] / ray.d[axis])
                step.append(1)
                out.append(self.n_voxels[axis])
            else:
                # handle ray with negative direction for voxel stepping
                next_crossing_t.append(
                    ray_t + (self._voxel_to_pos(pos[axis], axis) - \
                             grid_intersect[axis]) / ray.d[axis])
                delta_t.append(-self.width[axis] / ray.d[axis])
                step.append(-1)
                out.append(-1)

        # acquire a READ lock
        self.rw_lock.acquire_read()

        # walk ray through voxel grid
        hit_something = False
        while(True):
            # check for intersection in current voxel and advance to next
            voxel = self.voxels[self._offset(pos[0], pos[1], pos[2])]
            if voxel is not None:
                hit_something |= voxel.intersect(ray, intersection, self.rw_lock)

            # advance to next voxel

            # find step_axis for stepping to next voxel
            # don't use shift comparisons as it's slower than branching
            # in python (see /timing/minimum.py)
            if next_crossing_t[0] < next_crossing_t[1]:
                if next_crossing_t[0] < next_crossing_t[2]:
                    step_axis = 0
                else:
                    step_axis = 2
            else:
                if next_crossing_t[1] < next_crossing_t[2]:
                    step_axis = 1
                else:
                    step_axis = 2
            if ray.maxt < next_crossing_t[step_axis]:
                break
            pos[step_axis] += step[step_axis]
            if pos[step_axis] == out[step_axis]:
                break
            next_crossing_t[step_axis] += delta_t[step_axis]

        # release lock
        self.rw_lock.release()

        return hit_something

    def intersect_p(self, ray):
        """Return True if the ray intersects any primitive."""
        # check ray against overall grid bounds
        if self.bounds.inside(ray(ray.mint)):
            ray_t = ray.mint
        else:
            intersected, t0, t1 = self.bounds.intersect_p(ray)
            if not intersected:
                self.rw_lock.release()
                return False
            ray_t = t0

        grid_intersect = ray(ray_t)

        # set up 3D DDA (Digital Differential Analyzer) for ray
        pos = []
        next_crossing_t = []
        delta_t = []
        step = []
        out = []
        for axis in range(3):
            # compute current voxel for axis
            pos.append(self._pos_to_voxel(grid_intersect, axis))
            if ray.d[axis] == 0.0:
                next_crossing_t.append(float('inf'))
                delta_t.append(float('inf'))
                step.append(1)
                out.append(self.n_voxels[axis])
            elif ray.d[axis] > 0.0:
                # handle ray with positive direction for voxel stepping
                next_crossing_t.append(
                    ray_t + (self._voxel_to_pos(pos[axis]+1, axis) - \
                             grid_intersect[axis]) / ray.d[axis])
                delta_t.append(self.width[axis] / ray.d[axis])
                step.append(1)
                out.append(self.n_voxels[axis])
            else:
                # handle ray with negative direction for voxel stepping
                next_crossing_t.append(
                    ray_t + (self._voxel_to_pos(pos[axis], axis) - \
                             grid_intersect[axis]) / ray.d[axis])
                delta_t.append(-self.width[axis] / ray.d[axis])
                step.append(-1)
                out.append(-1)

        # acquire a READ lock
        self.rw_lock.acquire_read()

        # walk grid for shadow ray
        while(True):
            # check for intersection in current voxel and advance to next
            voxel = self.voxels[self._offset(pos[0], pos[1], pos[2])]
            if voxel and voxel.intersect_p(ray, self.rw_lock):
                self.rw_lock.release()
                return True

            # advance to next voxel

            # find step_axis for stepping to next voxel
            # don't use shift comparisons as it's slower than branching
            # in python (see /timing/minimum.py)
            if next_crossing_t[0] < next_crossing_t[1]:
                if next_crossing_t[0] < next_crossing_t[2]:
                    step_axis = 0
                else:
                    step_axis = 2
            else:
                if next_crossing_t[1] < next_crossing_t[2]:
                    step_axis = 1
                else:
                    step_axis = 2
            if ray.maxt < next_crossing_t[step_axis]:
                break
            pos[step_axis] += step[step_axis]
            if pos[step_axis] == out[step_axis]:
                break
            next_crossing_t[step_axis] += delta_t[step_axis]

        # release lock
        self.rw_lock.release()

        return False

    def _pos_to_voxel(self, point, axis):
        """Convert a 1D position into a voxel index."""
        v = int((point[axis] - self.bounds.p_min[axis]) *
                self.inv_width[axis])
        return clamp(v, 0, self.n_voxels[axis]-1)

    def _voxel_to_pos(self, index, axis):
        """Convert a voxel index to a 1D position."""
        return self.bounds.p_min[axis] + index * self.width[axis]

    def _offset(self, x, y, z):
        """Compute a voxel position based on its x, y, z indexes."""
        return z*self.n_voxels[0]*self.n_voxels[1] + y*self.n_voxels[0] + x

