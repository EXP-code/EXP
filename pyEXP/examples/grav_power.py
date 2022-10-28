#!/usr/bin/env python3

"""
Make amplitude/power plots using coefficient files
"""

import os, sys
import pyEXP
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
import sys

def help(phrase: str) -> None:
   """Print some usage info"""
   print(phrase)

def main() -> int:
    """Make the power plot"""

    if (len(sys.argv)==1):
        help(sys.argv[0] + ": you need at least one coefficient file on the command line")
        return 1

    plt.rcParams['figure.figsize'] = [12, 9]

    for n in range(1, len(sys.argv)):
        if not exists(sys.argv[n]): continue
        coefs = pyEXP.coefs.Coefs.factory(sys.argv[n])
        power = coefs.Power()
        for i in range(0, power.shape[1]):
            plt.semilogy(coefs.Times(), power[:,i], label='{}'.format(i))
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel(r'power=$[\sum_{n}|a_{mn}|^2]^{1/2}$')
        plt.title("Coefficient file="+sys.argv[n])
        plt.show()

    return 0

if __name__ == '__main__':
    sys.exit(main())
