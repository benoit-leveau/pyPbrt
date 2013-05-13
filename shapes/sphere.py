"""Sphere Class."""

import math

from core.pbrt import clamp, quadratic
from core.shape import Shape
from core.geometry import BBox, Point, Vector, Normal, Ray, \
     cross, dot, distance_squared, normalize
from core.diffgeom import DifferentialGeometry
from core.geometry import coordinate_system
# from core.monte_carlo import uniform_sample_sphere, uniform_sample_cone


class Sphere(Shape):
    
    """Class describing a Sphere."""

    def __init__(self, object_to_world, world_to_object, reverse_orientation,
                 radius, z0, z1, phi_max):
        """Default constructor for Sphere."""
        super(Sphere, self).__init__(object_to_world, world_to_object,
                                     reverse_orientation)
        self.radius = float(radius)
        self.z_min = clamp(min(z0, z1), -radius, radius)
        self.z_max = clamp(max(z0, z1), -radius, radius)
        self.theta_min = math.acos(clamp(self.z_min/radius, -1.0, 1.0))
        self.theta_max = math.acos(clamp(self.z_max/radius, -1.0, 1.0))
        self.phi_max = math.radians(clamp(phi_max, 0.0, 360.0))
        
    def object_bound(self):
        """Return bounding box in object space."""
        return BBox(Point(-self.radius, -self.radius, self.z_min),
                    Point( self.radius,  self.radius, self.z_max))

    def intersect(self, r):
        """Intersect the ray with the shape."""
        # Transform _Ray_ to object space
        ray = self.world_to_object(r)

        # Compute quadratic sphere coefficients
        A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z
        B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z)
        C = ray.o.x*ray.o.x + ray.o.y*ray.o.y + \
                  ray.o.z*ray.o.z - self.radius*self.radius

        # Solve quadratic equation for _t_ values
        found, t0, t1 = quadratic(A, B, C)
        if not found:
            return False, float('inf'), 0.0, None

        # Compute intersection distance along ray
        if (t0 > ray.maxt or t1 < ray.mint):
            return False, 0.0, 0.0, None
        t_hit = t0
        if (t0 < ray.mint):
            t_hit = t1
            if (t_hit > ray.maxt):
                return False, float('inf'), 0.0, None

        # Compute sphere hit position and $\phi$
        phi_t = ray(t_hit)
        if (phi_t.x == 0.0 and phi_t.y == 0.0):
            phi_t.x = 1e-5 * self.radius
        phi = math.atan2(phi_t.y, phi_t.x)
        if (phi < 0.0):
            phi += 2.0 * math.pi

        # Test sphere intersection against clipping parameters
        if ((self.z_min > -self.radius and phi_t.z < self.z_min) or \
            (self.z_max <  self.radius and phi_t.z > self.z_max) or \
            phi > self.phi_max):
            if (t_hit == t1):
                return False, float('inf'), 0.0, None
            if (t1 > ray.maxt):
                return False, float('inf'), 0.0, None
            t_hit = t1;
            # Compute sphere hit position and $\phi$
            phi_t = ray(t_hit);
            if (phi_t.x == 0.0 and phi_t.y == 0.0):
                phi_t.x = 1e-5 * self.radius;
            phi = math.atan2(phi_t.y, phi_t.x)
            if (phi < 0.0):
                phi += 2.0*math.pi
            if ((self.z_min > -self.radius and phi_t.z < self.z_min) or \
                (self.z_max <  self.radius and phi_t.z > self.z_max) or \
                phi > self.phi_max):
                return False, float('inf'), 0.0, None

        # Find parametric representation of sphere hit
        u = phi / self.phi_max
        theta = math.acos(clamp(phi_t.z / self.radius, -1.0, 1.0))
        v = (theta - self.theta_min) / (self.theta_max - self.theta_min)

        # Compute sphere $\dpdu$ and $\dpdv$
        zradius = math.sqrt(phi_t.x*phi_t.x + phi_t.y*phi_t.y)
        inv_z_radius = 1.0 / zradius
        cos_phi = phi_t.x * inv_z_radius
        sin_phi = phi_t.y * inv_z_radius
        dpdu = Vector(-self.phi_max * phi_t.y, self.phi_max * phi_t.x, 0)
        dpdv = (self.theta_max-self.theta_min) * \
               Vector(phi_t.z * cos_phi,
                      phi_t.z * sin_phi,
                      -self.radius * math.sin(theta))

        # Compute sphere $\dndu$ and $\dndv$
        d2Pduu = -self.phi_max * self.phi_max * Vector(phi_t.x, phi_t.y, 0)
        d2Pduv = (self.theta_max - self.theta_min) * phi_t.z * self.phi_max * \
                 Vector(-sin_phi, cos_phi, 0.0)
        d2Pdvv = -(self.theta_max - self.theta_min) * \
                 (self.theta_max - self.theta_min) * \
                 Vector(phi_t.x, phi_t.y, phi_t.z)

        # Compute coefficients for fundamental forms
        E = dot(dpdu, dpdu);
        F = dot(dpdu, dpdv);
        G = dot(dpdv, dpdv);
        N = normalize(cross(dpdu, dpdv))
        e = dot(N, d2Pduu)
        f = dot(N, d2Pduv)
        g = dot(N, d2Pdvv)

        # Compute $\dndu$ and $\dndv$ from fundamental form coefficients
        invEGF2 = 1.0 / (E*G - F*F)
        dndu = Normal.from_vector((f*F - e*G) * invEGF2 * dpdu + \
                                  (e*F - f*E) * invEGF2 * dpdv)
        dndv = Normal.from_vector((g*F - f*G) * invEGF2 * dpdu + \
                                  (f*F - g*E) * invEGF2 * dpdv)

        # Initialize _DifferentialGeometry_ from parametric information
        o2w = self.object_to_world
        dg = DifferentialGeometry.from_intersection(
            o2w(phi_t), o2w(dpdu), o2w(dpdv), o2w(dndu), o2w(dndv), u, v, self)

        # Compute _rayEpsilon_ for quadric intersection
        ray_epsilon = 5e-4 * t_hit
        return True, t_hit, ray_epsilon, dg

    def intersect_p(self, r):
        """Return True if ray intersects the shape."""
        # Transform _Ray_ to object space
        ray = self.world_to_object(r)

        # Compute quadratic sphere coefficients
        A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z
        B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z)
        C = ray.o.x*ray.o.x + ray.o.y*ray.o.y + \
            ray.o.z*ray.o.z - self.radius*self.radius

        # Solve quadratic equation for _t_ values
        found, t0, t1 = quadratic(A, B, C)
        if not found:
            return False

        # Compute intersection distance along ray
        if (t0 > ray.maxt or t1 < ray.mint):
            return False
        t_hit = t0
        if (t0 < ray.mint):
            t_hit = t1
            if (t_hit > ray.maxt):
                return False

        # Compute sphere hit position and $\phi$
        phi_t = ray(t_hit)
        if (phi_t.x == 0.0 and phi_t.y == 0.0):
            phi_t.x = 1e-5 * self.radius
        phi = math.atan2(phi_t.y, phi_t.x)
        if (phi < 0.0):
            phi += 2.0*math.pi

        # Test sphere intersection against clipping parameters
        if ((self.z_min > -self.radius and phi_t.z < self.z_min) or \
               (self.z_max <  self.radius and phi_t.z > self.z_max) or \
               phi > self.phi_max):
            if (t_hit == t1):
                return False
            if (t1 > ray.maxt):
                return False
            t_hit = t1
            # Compute sphere hit position and $\phi$
            phi_t = ray(t_hit)
            if (phi_t.x == 0.0 and phi_t.y == 0.0):
                phi_t.x = 1e-5 * self.radius
            phi = math.atan2(phi_t.y, phi_t.x)
            if (phi < 0.):
                phi += 2.0*math.pi
            if ((self.z_min > -self.radius and phi_t.z < self.z_min) or \
               (self.z_max <  self.radius and phi_t.z > self.z_max) or \
                phi > self.phi_max):
                return False
        return True

    def area(self):
        """Return the area of the shape."""
        return self.phi_max * self.radius * (self.z_max-self.z_min)

    def sample(self, u1, u2):
        """Sample the shape."""
        raise Exception("check_next_line")
        p = Point(0, 0, 0) + self.radius * 1.0 # uniform_sample_sphere(u1, u2)
        ns = normalize(self.object_to_world(Normal(p.x, p.y, p.z)))
        if (self.reverse_orientation):
            ns *= -1.0
        return self.object_to_world(p), ns

    def sample_p(self, p, u1, u2):
        """Sample at point p."""
        # Compute coordinate system for sphere sampling
        p_center = self.object_to_world(Point(0, 0, 0))
        wc = normalize(p_center - p)
        wc_x, wc_y = coordinate_system(wc)

        # Sample uniformly on sphere if $\pt{}$ is inside it
        if (distance_squared(p, p_center) - self.radius*self.radius) < 1e-4:
            return self.sample(u1, u2)

        # Sample sphere uniformly inside subtended cone
        sin_theta_max2 = self.radius*self.radius / distance_squared(p, p_center)
        cos_theta_max = math.sqrt(max(0.0, 1.0 - sin_theta_max2))
        raise Exception("next_line")
        # r = Ray(p, uniform_sample_cone(u1, u2, cos_theta_max, wcX, wcY, wc), 1e-3)
        r = Ray(p)
        intersect, t_hit, ray_epsilon, dg_sphere = self.intersect(r)
        if not intersect:
            t_hit = dot(p_center - p, normalize(r.d))
        ps = r(t_hit)
        ns = Normal(normalize(ps - p_center))
        if (self.reverse_orientation):
            ns *= -1.0
        return ps, ns

    def pdf_wi(self, p, wi):
        """Intersect sample ray with area light geometry."""
        p_center = self.object_to_world(Point(0, 0, 0))
        # Return uniform weight if point inside sphere
        if (distance_squared(p, p_center) - self.radius*self.radius) < 1e-4:
            return Shape.pdf_wi(self, p, wi)

        # Compute general sphere weight
        sin_theta_max2 = self.radius*self.radius / distance_squared(p, p_center)
        cos_theta_max = math.sqrt(max(0.0, 1.0 - sin_theta_max2))
        raise Exception("next_line")
        # return uniform_cone_pdf(cos_theta_max)
        return 0.0
    
    def __str__(self):
        """Return a string describing the sphere."""
        return "Sphere (%f, %f, %f, %f, %f, %f)" % (self.radius,
                                                    self.phi_max,
                                                    self.z_min,
                                                    self.z_max,
                                                    self.theta_min,
                                                    self.theta_max)
