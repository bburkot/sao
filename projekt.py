# -*- coding: utf-8 -*-
__author__ = 'Blazej'

import os
import sys
import matplotlib.pyplot as plt
import time
import numpy
import math

in_dir = "prime_numbers"
data = "data"

from_number = 2

to_number = 11

def simplify_data(in_dir, out_dir):
    i = 0
    for filename in os.listdir(in_dir):
        out = open(out_dir + "/" + filename, "w")
        for line in open(in_dir + "/" + filename):
            try:
                oLine = ""
                for number in line.split():
                    oLine = oLine + str(int(number)) + ";"

                if len(oLine) > 0:
                    out.write(oLine)
                    out.write("\n")
            except ValueError:
                pass
        out.close()
#simplify_data(in_dir, data)

hl = None

xdata = []
ydata = []

def create_plot(x, y):
    plt.hold(True)
    plt.ion()
    plt.figure()
    plt.plot(x, y, 'ro')
    plt.grid(True)
    plt.title('title')
    plt.ylabel('yTitle')
    plt.xlabel('xTitle')
    plt.axis([0, 7, 0, 6]) # [x0, x1, y0, y1]
    plt.savefig('result/out_' + time.strftime("%Y.%m.%d %H.%M.%S") + '.png')
    plt.draw()

def update_plot(x, y):
    plt.plot(x, y, 'ro')
    plt.draw()

def numberToVector(nr):

    pass

def work():
    x = [1,2,3,4]
    y = [1,2,3,4]
    create_plot(x, y)
    time.sleep( 1 )

    update_plot([5],[5])
    plt.ioff()
    #time.sleep(5)

def checkQuadrick(key):
    structs = {
        'right_top' : {
            'values':[1,3,13,31,57,91,133,183,241,307,381],
            'a': 4, 'b':-2, 'c':1
        },
        'right_mid':{
            'values':[1,2,11,28,53,86,127,176,233,298,371],
            'a': 4, 'b':-3, 'c':1
        },
        'right_bottom' : {
            'values':[1,9,25,49,81,121,169,225,289,361,441],
            'a': 4, 'b':4, 'c':1
        },
        'left_top' : {
            'values':[1,5,17,37,65,101,145,197,257,325,401],
            'a': 4, 'b':0, 'c':1
        },
        'left_mid' : {
            'values':[1,6,19,40,69,106,151,204,265,334,411],
            'a': 4, 'b':1, 'c':1
        },
        'left_bottom' : {
            'values':[1,7,21,43,73,111,157,211,273,343,421],
            'a': 4, 'b':2, 'c':1
        },
        'y_top' : {
            'values':[1,4,15,34,61,96,139,190],
            'a': 4, 'b':-1, 'c':1
        },
        'y_bottom' : {
            'values':[1,8,23,46,77,116,163,218,281,352,431],
            'a': 4, 'b':3, 'c':1
        },
    }

    y = structs[key]['values']
    a = structs[key]['a']
    b = structs[key]['b']
    c = structs[key]['c']

    for x in xrange(len(y)):
        if y[x] != a*(x**2) + b*x + c:
            print "rozne dla", x, "powinno byc", y[x],"a jest", a*(x**2) + b*x + c

def getXY(number):  #tested to 1.0 mln
    if number == 1:
        return 0,0

    n = math.ceil(math.sqrt((number-1)/4.0))
    right_mid = 4*n**2 - 3*n + 1
    top = 4*n**2 - n + 1

    if abs(right_mid - number) < n:
        return n, number - right_mid

    if abs(top - number) <= n:
        return  top - number, n

    n = math.floor(math.sqrt((number-1)/4.0))
    left_mid = 4*n**2 + n + 1
    bottom = 4*n**2 + 3* n + 1

    if abs(left_mid - number) <= n:
        return -n, left_mid - number

    if abs(bottom - number) <= n:
        return number - bottom, -n

    print "--------"
    print number, 'n ceil', math.ceil(math.sqrt((number-1)/4.0)), 'n floor',math.floor(math.sqrt((number-1)/4.0))
    print "bottom", bottom
    print "left_mid", left_mid
    print "top", top
    print "right_mid",right_mid
    return 0,0

#checkQuadrick('left_bottom')
# start = time.time() # in sec
# print "work time","%.5f" % (float(time.time()) - start),"[s]"

def generate_spiral(n):
    """ Return the first n coordinates in the Ulam spiral.  Note that
        (xs[i], ys[i]) is the location of the integer i+1 (i.e. 1 is located
        at (xs[0], ys[0])).
        """

    dirn = [0, 1] # The direction the spiral is traveling: [dX, dY]
    xs = [0, 1]
    ys = [0, 0]

    for i in range(2, n):
        # Compute new coordinates
        x = xs[-1] + dirn[0]
        y = ys[-1] + dirn[1]

        # Add to spiral
        xs.append(x)
        ys.append(y)

        # Change direction at the 45-degree lines (except SE line, which is
        # one unit past).  Draw it if you're trying to visualize...
        if x == y and x > 0:
            dirn = [-1, 0]
        elif -1 * x == y and x < 0:
            dirn = [0, -1]
        elif x == y and x < 0:
            dirn = [1, 0]
        elif x - 1 == -1 * y and x > 0:
            dirn = [0, 1]

    return xs, ys

def generate_spiral2(n):
    dirn = [0, 1] # The direction the spiral is traveling: [dX, dY]
    points = [[0,0],[1,0]]
    for i in xrange(2, n):
        # Compute new coordinates
        x = points[-1][0] + dirn[0]
        y = points[-1][1] + dirn[1]

        points.append([x, y])

        if x == y and x > 0:
            dirn = [-1, 0]
        elif -1 * x == y and x < 0:
            dirn = [0, -1]
        elif x == y and x < 0:
            dirn = [1, 0]
        elif x - 1 == -1 * y and x > 0:
            dirn = [0, 1]

    return points

def testGetXY(n):
    if n < 2:
        return;
    points = generate_spiral2(n)
    for i in xrange(n):
        r = getXY(i+1)
        if r[0] != points[i][0] or r[1] != points[i][1]:
            print 'number',i + 1,'get',r,'from spiral', points[i]

def testGetXYv2(n):
    dirn = [0, 1] # The direction the spiral is traveling: [dX, dY]
    lastX = 1
    lastY = 0

    i = 2
    while i < n:
        r = getXY(i)
        if r[0] != lastX or r[1] != lastY:
            print 'number',i,'get',r,'from spiral', [lastX, lastY]

        x = lastX + dirn[0]
        y = lastY + dirn[1]

        lastX = x
        lastY = y

        if x == y and x > 0:
            dirn = [-1, 0]
        elif -1 * x == y and x < 0:
            dirn = [0, -1]
        elif x == y and x < 0:
            dirn = [1, 0]
        elif x - 1 == -1 * y and x > 0:
            dirn = [0, 1]

        i += 1
    print getXY(i)

testGetXY(10**9)



# while True:
#     try:
#         nr = int(raw_input(">"))
#         ret = getXY(nr)
#         print "[X,Y]", ret
#     except ValueError:
#         print "pass"




