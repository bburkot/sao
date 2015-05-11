# -*- coding: utf-8 -*-
__author__ = 'Blazej Burkot'

import os
import sys
import time
import numpy
import math
import sqlite3
from sklearn.mixture import GMM
import numpy as np
import random
import decimal

ZERO = decimal.Decimal(0)
ONE = decimal.Decimal(1)
MINUSONE = decimal.Decimal(-1)
TWO = decimal.Decimal(2)
THREE = decimal.Decimal(3)
FOUR = decimal.Decimal(4)

class Conf:
    gmm_components = 4
    gmm_covariance_type = 'diag'               #   'spherical', 'tied', 'diag', 'full'
    gmm_n_iter = 100
    genetic_genPopulation = '_generatePopulation1'
    genetic_rate = "_rate1"
    genetic_crossover = "_crossover1"
    genetic_mutataion = "_mutation1"
    genetic_selection = "_selection1"
    work_getLine = 'getLine1'

decimal.getcontext().prec = 50

class DBconn(object):
    __db_name__ = "primes.db"

    _max_prime = 982451653   # it should be updated together with update in table primes
    _sizeOfTablePrimes = 50000000 # it too

    def __init__(self):
        self.conn = sqlite3.connect(DBconn.__db_name__)

    def getPrimes(self, fromNumber, toNumber=None):
        if toNumber == None:
            cursor = self.conn.execute("SELECT prime FROM primes WHERE id = :id", {'id' : fromNumber})
            num = cursor.fetchone()[0]
            cursor.close()
            return num
        else:
            cursor = self.conn.execute("SELECT prime FROM primes WHERE id >= :from and id <= :to", {'from' : fromNumber, "to" : toNumber})
            primes = cursor.fetchall()
            for i in xrange(len(primes)):
                primes[i] = primes[i][0]
            cursor.close()
            return primes

    def getRange(self, number):
        if number >= self._max_prime:
            return self._sizeOfTablePrimes - 1e5, self._sizeOfTablePrimes
        else:
            nbits = math.log(number, 2)
            fromNum = math.pow(2, nbits - 0.99)
            toNum = math.pow(2, nbits + 0.99)
            start = end = -1
            if fromNum <= 0:
                start = 0
            else:
                cursor = self.conn.execute("SELECT id FROM primes WHERE prime >= :fromNum and prime < :number", {"fromNum" : fromNum, 'number' : number})
                start = cursor.fetchone()[0]
                cursor.close()
            if toNum >= self._max_prime:
                end = self._sizeOfTablePrimes
            else:
                cursor = self.conn.execute("SELECT max(id) FROM primes WHERE prime > :number and prime <= :toNum", {"toNum" : toNum, 'number' : number})
                end = cursor.fetchone()[0]
                cursor.close()
            return (start, end)

    def getRange2(self, number):
        nbits = math.log(number, 2)
        if number > DBconn._max_prime:
            fromNum = DBconn._sizeOfTablePrimes - 1e5
        else:
            fromNum = math.pow(2, nbits - 0.99)
        toNum = math.pow(2, nbits + 0.99)
        cursor = self.conn.execute("SELECT prime FROM primes WHERE prime >= :fromNum and prime < :toNum", {"fromNum" : fromNum, 'toNum' : toNum})
        primes = []
        for row in cursor:
            primes.append(row[0])
        # primes = cursor.fetchall()
        # for i in xrange(len(primes)):
        #     primes[i] = int(primes[i][0])
        cursor.close()
        return primes

    def _getMaxPrime(self):
        cursor = self.conn.execute("SELECT max(prime) FROM primes")
        return cursor.fetchone()[0]
    def _getSizeOfTablePrimes(self):
        cursor = self.conn.execute("SELECT count(*) FROM primes")
        return cursor.fetchone()[0]
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()
    def __del__(self):
        self.conn.close()

class Line(object):
    def __init__(self, A = None, B = None, C = None):
        self.A = decimal.Decimal(A)
        self.B = decimal.Decimal(B)
        self.C = decimal.Decimal(C)
        if self.A == 0 and self.B == 0:
            raise Exception('A == 0 and B == 0')

    def calcDistance(self, x0, y0):
        return ((ONE * self.A * x0 + self.B * y0 + self.C) / (self.A ** 2 + self.B ** 2).sqrt()).copy_abs()

    def calcValue(self, x):
        return (self.A * x + self.C) / self.B

    def __str__(self):
        return "(" + str(self.A) + "*x + " + str(self.B) + "*y + " + str(self.C) + " = 0 )"
    def __repr__(self):
        return self.__str__()

class UlamSpiral(object):       #tested to 1e9
    @staticmethod
    def getXY(number):
        if number == 1:
            return ZERO, ZERO

        if number <= 0:
            raise NameError("Error: UlamSpiral.toXY(" + str(number) +") failed")

        n1 = (decimal.Decimal(number) - ONE) / FOUR #  n = math.ceil(math.sqrt((number-1)/4.0))
        n1 = n1.sqrt()

        n = decimal.Decimal(n1.to_integral_exact(decimal.ROUND_CEILING))

        right_mid = FOUR*n**2 - THREE*n + ONE
        top = FOUR*n**2 - n + ONE

        if abs(right_mid - number) < n:
            return n, number - right_mid

        if abs(top - number) <= n:
            return  top - number, n

        n = decimal.Decimal(n1.to_integral_exact(decimal.ROUND_FLOOR))

        left_mid = FOUR*n**2 + n + ONE
        bottom = FOUR*n**2 + THREE* n + ONE

        if abs(left_mid - number) <= n:
            return -n, left_mid - number

        if abs(bottom - number) <= n:
            return number - bottom, -n

        raise NameError("Error: UlamSpiral.toXY(" + str(number) +") failed")

    @staticmethod
    def getNumber(x = None, y = None, vector=None):
        if vector != None:
            x = vector.X
            y = vector.Y

        if x == y == ZERO:
            return ONE
        elif abs(y) >= abs(x):
            if y > ZERO:
                return FOUR*y*y - y + ONE - x
            else:
                return FOUR*y*y - THREE*y + ONE + x
        else:
            if x > ZERO:
                return FOUR*x*x - THREE*x + ONE + y
            else:
                return FOUR*x*x - x + ONE - y

    # @staticmethod
    # def getNumber(vec):
    #     return UlamSpiral.getNumber(vec.X, vec.Y)

class Vector(object):
    def __init__(self, number = None, x=None, y=None ):
        if number:
            self.X, self.Y = UlamSpiral.getXY(number)
        elif x  != None and y != None:
            self.X, self.Y = x, y
        else:
            raise AttributeError()

    def getLength(self):
        return (self.X**2 + self.Y**2).sqrt()

    def __str__(self):
        return "<" + str(self.X) + "," + str(self.Y) +">"
    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if type(other) == type(self):
            return Vector(x=(self.X + other.X), y=(self.Y + other.Y))

    def __div__(self, other):
        if isinstance( other, ( int, long ) ):
            return Vector(x=(self.X / other), y=(self.Y / other))

class Work(object):
    def __init__(self):
        self.getLine = self.__getattribute__('getLine1')

    def getDistances(self, rsaVec, primeVec):
        lineRsa = self.getLine(rsaVec.X, rsaVec.Y)
        orthogonalLine = Line(MINUSONE / lineRsa.A, MINUSONE, primeVec.Y + (ONE * primeVec.X) / lineRsa.A )

        # print lineRsa
        # print orthogonalLine

        d1 = orthogonalLine.calcDistance(rsaVec.X, rsaVec.Y) / rsaVec.getLength()
        d2 = lineRsa.calcDistance(primeVec.X, primeVec.Y) / rsaVec.getLength()

        if primeVec.Y >= lineRsa.calcValue(primeVec.X):
            return [d1, d2]
        else:
            return [d1, -d2]

    def getVector(self, rsaVec, d1, d2):
        pass


    def getLine1(self, x, y):
        if x > 0:
            return Line(1, -1, y - x)
        else:
            return Line(-1, -1, y + x)

def test():
    db = DBconn()
    prime1 = db.getPrimes(100)
    prime2 = db.getPrimes(115)
    rsa = prime1 * prime2

    rsaVec = Vector(rsa)
    prime1Vec = Vector(prime1)

    work = Work('getLine1')
    print work.getDistances(rsaVec, prime1Vec)


class GeneticAlgorithm(object):
    target_val = 0.0        # lower better but must by >= 0
    population_size = int(1e4)
    probability_crossover = 0.7
    probability_mutation = 0.4
    add_new_population_percent = 0.001

    def __init__(self, numberToFactorize):
        #self.gmm = gmm
        self.number =  decimal.Decimal(numberToFactorize)

        x, y = UlamSpiral.getXY(self.number.sqrt().to_integral_exact(decimal.ROUND_CEILING))
        self.searchBorder = max(abs(x), abs(y)) + 1 # int(sqrt(number to factorize)) is inside of the border

        if (Conf.genetic_genPopulation == 'generatePopulation2'):
            self.gmm = self._generateGauseMixture(numberToFactorize)
        self.generatePopulation = self.__getattribute__(Conf.genetic_genPopulation)
        self.rate = self.__getattribute__(Conf.genetic_rate)
        self.crossover = self.__getattribute__(Conf.genetic_crossover)
        self.mutation = self.__getattribute__(Conf.genetic_mutataion)
        self.selection = self.__getattribute__(Conf.genetic_selection)

        self.fitness = sys.maxint
        self.results = []


    """ GENERATE POPULATION  """
    def _generatePopulation1(self, size):
        vectors = []
        for _ in xrange(size):
            x1, y1  = numpy.random.randint(-1 * self.searchBorder, self.searchBorder + 1, size=2)
            vectors.append(Vector(x=x1, y=y1))
        return vectors

    def _generatePopulation2(self, size):
        pass
    def _generateGauseMixture(number):
        print time.ctime(), "start generate gause mixture"
        db = DBconn()
        work = Work('getLine1')

        start, end = db.getRange(number)

        samples = np.random.randint(start, end, size=(Conf.generate_samples, 2))

        testData = []
        for el in samples:
            p1 = db.getPrimes(el[0])
            p2 = db.getPrimes(el[1])
            try:
                rsaVec = Vector(p1 * p2)
                prime1Vec = Vector(p1)
                prime2Vec = Vector(p2)
            except Exception as e:
                print e
                print p1, p2, p1 * p2, type(p1), type(p2),type(64888757)

                p1 = p1.item()
                p2 = p2.item()

                rsaVec = Vector(p1 * p2)
                prime1Vec = Vector(p1)
                prime2Vec = Vector(p2)

            testData.append(work.getDistances(rsaVec, prime1Vec))
            testData.append(work.getDistances(rsaVec, prime2Vec))

        gmm = GMM(n_components = Conf.gmm_components, covariance_type = Conf.gmm_covariance_type, n_iter = Conf.gmm_n_iter)
        # print gmm.fit(testData)
        #
        # print gmm.get_params()
        # print 'means'
        # print gmm.means_
        # print 'covaras'
        # print gmm.covars_
        # print 'weights'
        # print gmm.weights_
        # print 'sample'
        # print gmm.sample(10)

        db.conn.close()
        print time.ctime(), "end generate gause mixture"
        return gmm

    """ RATE """
    def _rate1(self, vec):
        number = UlamSpiral.getNumber(vector=vec)
        if number == self.number or number == 1:
            return sys.maxint

        value = self.number % number
        if value == 0:
            self.fitness = 0
            self.results.append(number)

        return value

    """ CROSSOVER """
    def _crossover1(self, vec1, vec2):
        shift = 5
        shared = 2**shift -1

        if random.random() > 0.5:
            x = ((vec1.X >> shift) << shift) + (vec2.X & shared)
        else:
            x = ((vec2.X >> shift) << shift) + (vec1.X & shared)

        if random.random() > 0.5:
            y = ((vec1.Y >> shift) << shift) + (vec2.Y & shared)
        else:
            y = ((vec2.Y >> shift) << shift) + (vec1.Y & shared)

        return (vec1 + vec2) /2 , Vector(x=x, y=y)

    """ MUTATION """
    def _mutation1(self, vec):
        x = int(vec.X)
        xbits = math.floor(numpy.log2(abs(x) + 1) + 1)
        change =  numpy.random.randint(xbits, size=numpy.random.randint(xbits))
        for nbit in change:
            x ^= 2 ** nbit

        y = int(vec.X)
        ybits = math.floor(numpy.log2(abs(y) + 1) + 1)
        change =  numpy.random.randint(ybits, size=numpy.random.randint(ybits))
        for nbit in change:
            y ^= 2 ** nbit
        return Vector(x=x, y=y)

    """ SELECTION """
    def _selection1(self, vectors):
        temp = []
        newPopulation = self._generatePopulation1(int(math.ceil(GeneticAlgorithm.population_size * GeneticAlgorithm.add_new_population_percent)))
        vectors += newPopulation
        for vec in vectors:
            temp.append([vec, self.rate(vec)])

        sortedIter = sorted(temp, key = lambda x: x[1])

        i = 0
        ret = []
        for item in sortedIter:
            if item[0].X <= self.searchBorder and item[0].Y <= self.searchBorder:
                ret.append(item[0])
                i += 1
                if i >= GeneticAlgorithm.population_size:
                    break
        if len(ret) < GeneticAlgorithm.population_size:
            newPopulation = self._generatePopulation1(math.ceil(GeneticAlgorithm.population_size - i))
            ret += newPopulation
        return ret
    """ END of prototypes of function """

    def run(self):
        print time.ctime(), "- start genetic alg"
        vectors = self.generatePopulation(GeneticAlgorithm.population_size)
        while True:
            vectors = self.selection(vectors)
            if GeneticAlgorithm.target_val >= self.fitness : break

            size = len(vectors)
            for i in xrange(size):
                if random.random() < GeneticAlgorithm.probability_crossover:
                    ch1, ch2 = self.crossover(vectors[i], vectors[ (i+1)%size])
                    vectors.append(ch1)
                    vectors.append(ch2)
                elif random.random() < GeneticAlgorithm.probability_mutation:
                    vectors.append(self.mutation(vectors[i]))
        print time.ctime(), "- end genetic alg"
        return self.fitness, self.results




# print time.ctime(), "start "
# generateGauseMixture(3846564564645645646546468798909808908908090890890768566354424214264576890000000000000000000000000000000000000000000000000000000000242352)
# print time.ctime(), "end "


alg = GeneticAlgorithm(5018719) #Factorization: 1823 * 2753 = 5018719
#alg = GeneticAlgorithm(695903367050368781) # 756067723 * 920424647 = 695903367050368781
alg.run()
#print alg.fitness
print alg.results




