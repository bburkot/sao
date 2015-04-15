# -*- coding: utf-8 -*-
__author__ = 'Blazej'

import os
import sys
import matplotlib.pyplot as plt
import time
import numpy
import math
import sqlite3

db_name = "primes.db"


def createDB():
    conn = sqlite3.connect(db_name)
    conn.execute("CREATE TABLE PRIMES (ID INTEGER PRIMARY KEY AUTOINCREMENT, PRIME INTEGER NOT NULL);")
    for nr in xrange(1,51):
        file = open("data/primes%d.txt" % nr);
        for line in file:
            for nr in line.split(";"):
                if nr != '\n':
                    conn.execute("INSERT INTO PRIMES (PRIME) VALUES (" + nr + ")");
        conn.commit()

    cursor = conn.execute("SELECT count(*) FROM primes ")
    print cursor.fetchone()[0]

    conn.close()

def getPrimes(fromNumber, toNumber=None):
    global conn
    if toNumber == None:
        cursor = conn.execute("SELECT prime FROM primes WHERE id = :id", {'id' : fromNumber})
        return cursor.fetchone()[0]
    else:
        primes = []
        cursor = conn.execute("SELECT prime FROM primes WHERE id >= :from and id <= :to", {'from' : fromNumber, "to" : toNumber})
        for row in cursor:
            primes.append(row[0])
        return primes


class Line:
    def __init__(self, A = None, B = None, C = None):
        self.A = A
        self.B = B
        self.C = C
        if self.A == 0 and self.B == 0:
            raise Exception('A == 0 and B == 0')
    def __str__(self):
        return str(self.A) + "*x + " + str(self.B) + "*y + " + str(self.C) + " = 0"
    def __repr__(self):
        return self.__str__()

    def calcDistance(self, x0, y0):
        return abs(1.0 * self.A * x0 + self.B * y0 + self.C) / math.sqrt(self.A ** 2 + self.B ** 2)

    def calcValue(self, x):
        return (self.A * x + self.C) / self.B

class Vector:
    def __init__(self, number):
        self.x1, self.y1 = Vector.__numberToVector__(number)

    def __str__(self):
        return "[" +str(self.x1) + ", " + str(self.y1) + "]"
    def __repr__(self):
        return self.__str__()

    def __countMatchingLine__(self):
        if self.x1 > 0:
            self.line = Line(1, -1, self.y1 - self.x1)
        else:
            self.line = Line(-1, -1, self.y1 + self.x1)

    @staticmethod
    def __numberToVector__(number):  #tested to 1.0 mln
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

    @staticmethod
    def createMainVector(number):
        v = Vector(number)
        v.__countMatchingLine__()
        return v

    def getVectors(self, vec):
        return self.getVectors2(vec.x1, vec.y1)

    def getVectors2(self, x0, y0):
        distBeetweenPointAndLine = self.line.calcDistance(x0,y0)
        distBetweenPoints = math.sqrt((x0 - self.x1)**2 + (y0 - self.y1)**2)
        dist = math.sqrt(distBetweenPoints**2 - distBeetweenPointAndLine**2)
        if y0 >= self.line.calcValue(x0):
            return [dist, distBeetweenPointAndLine]
        else:
            return [dist, -distBeetweenPointAndLine]

    def getLength(self):
        return math.sqrt(self.x1**2 + self.y1**2)

class Presentation:
    def __init__(self):
        self.dir = "result/start_time_" + time.strftime("%Y.%m.%d %H.%M.%S") + "/"
        os.makedirs(self.dir)
        self.plotNumber = 1

    def create_plot(self, x, y, title, xTitle, yTitle):
        plt.hold(True)
        plt.ion()
        figure = plt.figure()
        plt.scatter(x, y, s=0.1)
        plt.grid(True)
        plt.title(title)
        plt.ylabel(yTitle)
        plt.xlabel(xTitle)
       # plt.axis([-1, 2, -2, 2]) # [x0, x1, y0, y1]
        plt.savefig(self.dir + '/out' +str(self.plotNumber) + '_' + time.strftime("%Y.%m.%d %H.%M.%S") + '.png')
        self.plotNumber += 1
        #plt.draw()
        return figure

    def update_plot(self,x, y):
        plt.plot(x, y, 'ro')
        plt.savefig(self.dir + '/out_' + time.strftime("%Y.%m.%d %H.%M.%S") + '.png')
        plt.draw()

def gen_stats(count):

    dataLenOfTheRSANumber = []
    dataDistFirstNumber = []
    dataDistSecondNumber = []

    primes = getPrimes(1, count)

    for i in xrange(len(primes)):
        vecPrime1 = Vector(primes[i])
        for j in xrange(i+1, len(primes)):
            vecPrime2 = Vector(primes[j])
            rsa = primes[i] * primes[j]
            if numpy.absolute(numpy.log2(primes[i]) - numpy.log2(primes[j])) < 1:
                vecRsa = Vector.createMainVector(rsa)
                vesRsaPrime = Vector.createMainVector(math.sqrt(rsa))

                #RSA number
                #dataLenOfTheRSANumber.append(numpy.log2(rsa))
                dataLenOfTheRSANumber.append(numpy.log2(rsa))

                #prime 1 todo
                distOnLine, distFromLine = vesRsaPrime.getVectors(vecPrime1)
                dataDistFirstNumber.append(distFromLine)

                #prime 2 todo
                distOnLine, distFromLine = vesRsaPrime.getVectors(vecPrime2)
                dataDistSecondNumber.append(distFromLine)

    p = Presentation()
    p.create_plot(dataLenOfTheRSANumber, dataDistFirstNumber, 'title', 'Len of RSA number in bits', 'Distance from line first number')
    p.create_plot(dataLenOfTheRSANumber, dataDistSecondNumber, 'title', 'Len of RSA number in bits', 'Distance from line second number')

conn = sqlite3.connect(db_name)

#gen_stats(2000) # PS zwróć uwagę że są jakieś wzorki na wykresach

vecRsa = Vector.createMainVector(188198812920607963838697239461650439807163563379417382700763356422988859715234665485319060606504743045317388011303396716199692321205734031879550656996221305168759307650257059)
distOnLine, distFromLine = vecRsa.getVectors(Vector(472772146107435302536223071973048224632914695302097116459852171130520711256363590397527))
print distFromLine;


conn.close()


#pozostaje do sprawdzenia jeszcze moja metoda
#wpisac na sztywno liczby RSA i zobaczyc jak blisko bedzie ktorys z dzielnikow