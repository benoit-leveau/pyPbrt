#!/usr/bin/python

import sys

from core.geometry import Point, Vector
from core.transform import translate, look_at, Transform, rotate, AnimatedTransform
from core.primitive import GeometricPrimitive, TransformedPrimitive
from core.scene import Scene
from core.integrator import Integrator
from shapes.sphere import Sphere
from accelerators.grid import GridAccel
from renderers.sampler_renderer import SamplerRenderer
from integrators.whitted import WhittedIntegrator
from film.image import ImageFilm
from cameras.perspective import PerspectiveCamera
from samplers.random import RandomSampler
from filters.box import BoxFilter


def create_pyramid():
    primitives = []
    # create the sphere primitive
    material = None
    object_to_world = Transform()
    world_to_object = Transform()
    sphere = Sphere(object_to_world, world_to_object, False, 1.0, -1.0,
                    1.0, 360)
    primitive = GeometricPrimitive(sphere, material, None)
    
    # now create the instances
    # create the objects, with this layout (shown in 2D)
    #   O O O O  level 0, width 4, pos_x -6
    #    O O O   level 1, width 3, pos_x -4
    #     O O    level 2, width 2, pos_x -2
    #      O     level 3, width 1, pos_x 0

    for level in range(4):
        width_array = 4 - level
        start_pos = Point(-2.0*(3-level), -2.0*(3-level), level*2.0)
        for i in range(width_array):
            start_pos_i = start_pos + i * Vector(2.0, 0.0, 0.0)
            for j in range(width_array):
                pos = start_pos_i + j * Vector(0.0, 2.0, 0.0)
                transform = translate(pos).inverse()
                world_to_primitive = AnimatedTransform(transform, 0.0, transform, 1.0)
                instance = TransformedPrimitive(primitive, world_to_primitive)
                primitives.append(instance)
    return primitives


def create_simple_sphere():
    primitives = []
    material = None
    object_to_world = translate(Point(0,1,0))
    world_to_object = object_to_world.inverse()
    sphere = Sphere(object_to_world, world_to_object, False,
                    1.0, -1.0, 1.0, 360)
    primitive = GeometricPrimitive(sphere, material, None)
    primitives.append(primitive)
    return primitives

    
def create_scene(primitives):
    """Create a default scene."""

    # create the accelerator
    aggregate = GridAccel(primitives, False)

    # create the lights
    lights = []
    # light = PointLight()
    # lights.append(light)

    # create the scene
    scene = Scene(aggregate, lights, volume_region=None)

    return scene


def create_film(filename, width, height):
    image_filter = BoxFilter(0.5, 0.5)
    crop = [0.0, 1.0, 0.0, 1.0]
    open_window = False
    film = ImageFilm(width, height, image_filter, crop, filename, open_window)
    return film


def create_camera(film, cam_transform):
    """Create a perspective camera."""
    cam2world = cam_transform.inverse()
    frame = float(film.x_resolution) / float(film.y_resolution)
    if frame > 1.0:
        screen_window = [-frame, frame, -1.0, 1.0]
    else:
        screen_window = [-1.0, 1.0, -1.0/frame, 1.0/frame]
    s_open = 0.0
    s_close = 0.5
    lens_radius = 0.0
    focal_distance = 10.0
    fov = 22.0

    camera = PerspectiveCamera(cam2world,
                               screen_window,
                               s_open, s_close,
                               lens_radius, focal_distance,
                               fov, film)

    return camera


def create_renderer(camera):
    """Create a renderer."""
    # surface integrator
    whitted_integrator = WhittedIntegrator(3)

    # volume integrator
    volume_integrator = Integrator()

    # sampler
    sampler = RandomSampler(0, camera.film.x_resolution-1,
                            0, camera.film.y_resolution-1,
                            16,
                            0.0, 0.5
                            )

    # create the renderer
    sampler_renderer = SamplerRenderer(sampler,
                                       camera,
                                       whitted_integrator,
                                       volume_integrator
                                       )

    return sampler_renderer


def main():
    """Main function."""

    if len(sys.argv) < 5:
        print "Error. Specify filename, sene, width, height."
        return -1
    filename = sys.argv[1]
    scene = sys.argv[2]
    width = int(sys.argv[3])
    height = int(sys.argv[4])

    if scene == "sp":
        primitives = create_simple_sphere()
        cam_pos = Point(0, 3, 8)
        cam_look = Point(0, 0.8, 0)
        cam_up = Vector(0, 1, 0)
    elif scene == "py":
        primitives = create_pyramid()
        cam_pos = Point(10, 10, 4)
        cam_look = Point(0, 0, 2)
        cam_up = Vector(0, 0, 1)
    else:
        print "unknown scene"
        return -1
        
    cam_transform = look_at(cam_pos, cam_look, cam_up)
    scene = create_scene(primitives)
    film = create_film(filename, width, height)
    camera = create_camera(film, cam_transform)
    renderer = create_renderer(camera)

    # launch render
    renderer.render(scene)

    return 0

if __name__ == '__main__':
    sys.exit(main())
