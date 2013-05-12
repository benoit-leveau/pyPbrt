#!/usr/bin/python

import random
import time

nb_calls = 10000
size = 100

cmp = [2,1,2,1,2,2,0,0]

def shift_min(a, b, c):
    return cmp[((a<b)<<2) + ((a<c)<<1) + (b<c)]

def branch_min(a, b, c):
    if a<b:
        if a<c:
            return 0
        else:
            return 2
    if b<c:
        return 1
    return 2

def run_test(is_float=False):

    values = []
    for i in xrange(size):
        if is_float:
            a = random.random() * 20.0
            b = random.random() * 20.0
            c = random.random() * 20.0
        else:
            a = random.randint(0,20)
            b = random.randint(0,20)
            c = random.randint(0,20)
        values.append((a,b,c))
    
    t1 = time.time()

    for i in xrange(nb_calls):
        a, b, c = values[i%size]
        shift_min(a, b, c)

    t2 = time.time()

    for i in xrange(nb_calls):
        a, b, c = values[i%size]
        branch_min(a, b, c)

    t3 = time.time()

    return t1, t2, t3

def test():
    t1, t2, t3 = run_test()
    print "= INTEGER ="
    print "%d calls" % (nb_calls)
    print "shift_comp   %.2fns" % ((t2-t1)/float(nb_calls)*1e6)
    print "branch_comp  %.2fns" % ((t3-t2)/float(nb_calls)*1e6)

    t1, t2, t3 = run_test(True)
    print "= FLOAT ="
    print "%d calls" % (nb_calls)
    print "shift_comp   %.2fns" % ((t2-t1)/float(nb_calls)*1e6)
    print "branch_comp  %.2fns" % ((t3-t2)/float(nb_calls)*1e6)
    

if __name__ == '__main__':
    test()

