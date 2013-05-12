#!/usr/bin/python

import sys

from core.geometry import Point, Vector
from core.transform import translate, look_at
from core.primitive import GeometricPrimitive
from core.scene import Scene
from core.integrator import Integrator
from shapes.sphere import Sphere
from accelerators.grid import GridAccel
from renderers.sampler_renderer import SamplerRenderer
from integrators.whitted import WhittedIntegrator
from film.image import ImageFilm
from cameras.perspective import PerspectiveCamera


def create_scene():
    """Create a default scene."""
    
    # create the objects, with this layout (shown in 2D)
    #   O O O O  level 0, width 4, pos_x -6
    #    O O O   level 1, width 3, pos_x -4
    #     O O    level 2, width 2, pos_x -2
    #      O     level 3, width 1, pos_x 0
    primitives = []
    material = None
    for level in range(4):
        width_array = 4 - level
        start_pos = Point(-2.0*(3-level), -2.0*(3-level), level*2.0)
        for i in range(width_array):
            start_pos_i = start_pos + i * Vector(2.0, 0.0, 0.0)
            for j in range(width_array):
                pos = start_pos_i + j * Vector(0.0, 2.0, 0.0)
                object_to_world = translate(pos)
                world_to_object = object_to_world.inverse()
                sphere = Sphere(object_to_world, world_to_object,
                                False, 1.0, -1.0, 1.0, 360)
                primitive = GeometricPrimitive(sphere, material, None)
                primitives.append(primitive)

    # create the accelerator
    grid = GridAccel(primitives, False)

    # create the lights
    lights = []
    # light = PointLight()
    # lights.append(light)

    # create the scene
    scene = Scene(grid, lights, volume_region=None)

    return scene


def create_film(filename):
    width = 640
    height = 480
    image_filter = None
    crop = [0,0,0,0]
    open_window = False
    film = ImageFilm(width, height, image_filter, crop, filename, open_window)
    return film


def create_camera(film):
    """Create a perspective camera."""
    cam_pos = Point(20, 20, 5)
    cam_look = Point(0, 0, 0)
    cam_up = Vector(0, 0, 1) 
    cam_to_world = look_at(cam_pos, cam_look, cam_up)
    screen_window = [0, 639, 0, 479]
    s_open = 0.0
    s_close = 0.5
    lens_radius = 35.0
    focal_distance = 35.0
    fov = 40.0
    camera = PerspectiveCamera(cam_to_world,
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
    sampler = None

    # create the renderer
    sampler_renderer = SamplerRenderer(sampler,
                                       camera,
                                       whitted_integrator,
                                       volume_integrator
                                       )

    return sampler_renderer


def main():
    """Main function."""

    if len(sys.argv) < 2:
        print "Error. Specify filename."
        return -1
    filename = sys.argv[1]

    scene = create_scene()
    film = create_film(filename)
    camera = create_camera(film)
    renderer = create_renderer(camera)

    # launch render
    renderer.render(scene)

    return 0

if __name__ == '__main__':
    sys.exit(main())
