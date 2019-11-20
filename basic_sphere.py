#!/usr/env/python
#//taken from https://9bitscience.blogspot.com/2013/07/raymarching-distance-fields_14.html
import numpy as np

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
       return v
    return v / norm

resolution = [80,50]
maxsteps = 32
max_dist = 10
gl_FragCoord = [300,300]

# main
right = [1, 0, 0]
up = [0,1,0]
#print rd
color = 0
total_steps = 0

# raymarch a sphere
for xx in range(resolution[1]):
    for yy in range(resolution[0]):
        u = xx * 2.0 / resolution[1] - 1.0
        v = yy * 2.0 / resolution[0] - 1.0
        ro = np.multiply(right,u) + np.multiply(up,v*1.77)
        rd = normalize(np.cross(right, up))
        color = 0
        t = 0
        for i in range(maxsteps):
            total_steps += 1
            p = ro + np.multiply(rd, t);
            #print p
            d = np.linalg.norm(p) - 0.5
            if (d > max_dist):
                break
            #print d
            if (d < 0.01):
                #print "white !"
                color = 1
                break
            t += d
        print color,
    print ""

print "total_steps: %d" % total_steps
