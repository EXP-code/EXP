R"(#!/usr/bin/python3

# Try to prevent ChiantiPy from grabbing the mouse                              
#
import matplotlib as mpl
mpl.use('Agg')

# Now import ChiantiPy
#
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import math
from optparse import OptionParser

elementList = {1:"h", 2:"he", 3:"li", 4:"be", 5:"b", 6:"c", 7:"n", 8:"o",9:"f",10:"ne", 11:"na", 12:"mg", 13:"al", 14:"si", 15:"p", 16:"s", 17:"cl", 18:"ar", 19:"k", 20:"ca", 21:"sc", 22:"ti", 23:"v", 24:"cr", 25:"mn", 26:"fe", 27:"co", 28:"ni", 29:"cu", 30:"zn", 31:"ga", 32:"ge", 33:"as", 34:"se", 35:"br", 36:"kr"}

def getRecomb(ofile='recomb.dat', Z=1, tmin=1.0e+03, tmax=1.0e+07, num=200):
    """ Use ChiantiPy to get the recombination rates"""

    if Z not in elementList:
        print("Element {} is not in list".format(Z))
        return

    file = open(ofile, 'w')

    Tmin = math.log(tmin)
    Tmax = math.log(tmax)
    dT   = (Tmax - Tmin)/(num-1)
    temp = []
    for i in range(num+1): temp.append(math.exp(Tmin + dT*i))

    data = []

    for z in range(2,Z+2):
        name = '{}_{}'.format(elementList[Z], z)
        v = ch.ion(name, temperature=temp)
        v.recombRate()
        data.append(v.RecombRate['rate'])

    for i in range(num+1):
        print("{:16.8e}".format(temp[i]), end='', file=file)
        for z in range(2,Z+2):
            print("{:16.8e}".format(data[z-2][i]), end='', file=file)
        print(file=file)

def main():
        """
        Get the ionization equilibrium for some number of elements
        Print to a file
        """
        usage  = "usage: %prog [options] file"
        parser = OptionParser(usage)

        parser.add_option("-Z", "--element", default=1,
                      action="store", type="int", dest="Z",
                      help="atomic element")
        parser.add_option("-n", "--number", default=100,
                      action="store", type="int", dest="number",
                      help="number of temperatures")
        parser.add_option("-t", "--minT", default=1000.0,
                      action="store", type="float", dest="minT",
                      help="minimum temperature")
        parser.add_option("-T", "--maxT", default=10000000.0,
                      action="store", type="float", dest="maxT",
                      help="maximum temperature")
        parser.add_option("-o", "--output", default='element.dat',
                      action="store", type="string", dest="ofile",
                      help="output file")

        (opt, args) = parser.parse_args()

        # Call the main routine
        getRecomb(opt.ofile, opt.Z, tmin=opt.minT, tmax=opt.maxT, num=opt.number)

if __name__ == "__main__":
        main()

)"
