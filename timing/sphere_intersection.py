#!/usr/bin/python

import random
import time

from core.transform import Transform, translate, scale
from shapes.sphere import Sphere
from core.geometry import Ray, Vector, Point

size = 100
nb_calls = 10000

def time_sphere_intersection():
    # create a transform
    o2w = translate(Vector(10,0,0)) * scale(1.3, 1.8, 2.0)
    w2o = o2w.inverse()

    # create the sphere
    sphere = Sphere(o2w, w2o, False, 1.0, -1.0, 1.0, 360)

    # create a large amount of rays,
    # choose so that half of them will intersect the ray

    positions = [Point(random.randint(0,100),
                       random.randint(0,100),
                       random.randint(0,100)
                       ) for i in range(size)]

    ray = Ray(Point(0,0,0), Vector(1.0, 1.0, 1.0))
    vectors = []
    for i in xrange(size):
        position = positions[i]
        if i%2 == 0:
            # make sure this ray hit the sphere
            vector = sphere.object_to_world(Point(0, 0, 0)) - position
            vector /= float(random.randint(1,10))
        else:
            # construct a random vector
            vector = Vector((random.random()-0.5)*random.randint(1,5),
                            (random.random()-0.5)*random.randint(1,5),
                            (random.random()-0.5*random.randint(1,5)))
                            
        vectors.append(vector)

    intersections = 0
    t1 = time.time()
    for i in xrange(nb_calls):
        ray.o = positions[i%size]
        ray.d = vectors[i%size]
        if sphere.intersect_p(ray):
            intersections += 1

    t2 = time.time()
    for i in xrange(nb_calls):
        ray.o = positions[i%size]
        ray.d = vectors[i%size]
        sphere.intersect(ray)

    t3 = time.time()

    print "%d calls, %d intersections" % (nb_calls, intersections)
    print "Sphere.intersect_p() %.2fms" % ((t2-t1)/float(nb_calls)*1000.0)
    print "Sphere.intersect()   %.2fms" % ((t3-t2)/float(nb_calls)*1000.0)

if __name__ == '__main__':
    time_sphere_intersection()

