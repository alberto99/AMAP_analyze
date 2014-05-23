#!/usr/bin/env python

#import graphlib
from collections import Counter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import numpy
import random
import sys
# Force matplotlib to not use any Xwindows backend.
import matplotlib.cm as cm


def add_noise(data, width):
    n = len(data)
    for i in range(n):
        data[i] += random.uniform(-width, width)
    return data


def get_data_for_target(target):
    rmsd = []
    rmsdt = []
    volume = []
    results = target.results
    for p in results:
        rmsd.append(p.rmsd_native)
        rmsdt.append(p.rmsd_template)
        volume.append(p.volume)
    return numpy.array(rmsd),numpy.array(rmsdt),numpy.array(volume)

def return_consensus(data):
    #   AlphaHelix { set assign "H"}
    #   310Helix { set assign "G"}
    #   PiHelix { set assign "I"}
    #   Strand { set assign "E"}
    #   Bridge { set assign "B"}
    #   Coil { set assign "C"}
    #   Turn* { set assign "T"}
    #   Gamma* { set assign "T"}
    dictio = {"H":0,"G":1,"I":2,"E":3,"B":4,"C":5,"T":6}
    hel = {0:"H",1:"G",2:"I",3:"E",4:"B",5:"C",6:"T"}

    sse = []
    ahelix=[]
    helix310=[]
    pihelix=[]
    strand=[]
    bridge=[]
    coil=[]
    turn=[]

    consensus = []
    for i in range(len(data[0,:])):
        a = list(data[:,i])
        num = []
        num.append(a.count("H"))
        num.append(a.count("G"))
        num.append(a.count("I"))
        num.append(a.count("E"))
        num.append(a.count("B"))
        num.append(a.count("C"))
        num.append(a.count("T"))
        index = numpy.argmax(numpy.array(num))
        consensus.append(hel[index])
    return consensus


def main():

    data = numpy.genfromtxt(sys.argv[1],dtype='str')
    consensus = return_consensus(data)
    print "".join(consensus)

if __name__ == '__main__':
    main()

