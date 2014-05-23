#!/usr/bin/env python

#import graphlib
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

def read_file(filen="",l=[2,7,12]):
    """
    read ptraj output file: time, rmsd
    """
    fn = open(filen,"r")
    #   AlphaHelix { set assign "H"}
    #   310Helix { set assign "G"}
    #   PiHelix { set assign "I"}
    #   Strand { set assign "E"}
    #   Bridge { set assign "B"}
    #   Coil { set assign "C"}
    #   Turn* { set assign "T"}
    #   Gamma* { set assign "T"}
    dictio = {"H":0,"G":1,"I":2,"E":3,"B":4,"C":5,"T":6}

    sse = []
    ahelix=[]
    helix310=[]
    pihelix=[]
    strand=[]
    bridge=[]
    coil=[]
    turn=[]

    count = 0
    for line in fn:
        if count < 999:
            count += 1
            continue
        if count > 10000:
            break
        count += 1
        a = line.split()
        if len(l) < 1:
            l = range(0,len(a))
        for i in l:
            ahelix.append(a[i].count("H"))
            helix310.append(a[i].count("G"))
            pihelix.append(a[i].count("I"))
            strand.append(a[i].count("E"))
            bridge.append(a[i].count("B"))
            coil.append(a[i].count("C"))
            turn.append(a[i].count("T"))

        elem = []
        for i in a:
            elem.append(dictio[i])      
        sse.append(numpy.array(elem))

    a_hel = numpy.average(ahelix)*100
    a_310 = numpy.average(helix310)*100
    a_pi = numpy.average(pihelix)*100
    a_beta = numpy.average(strand)*100
    a_bridge = numpy.average(bridge)*100
    a_coil = numpy.average(coil)*100
    a_turn = numpy.average(turn)*100
    #print "name a_hel,a_310,a_pi,a_beta,a_bridge,a_coil,a_turn"
    print "%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f" % (a_hel,a_310,a_pi,a_beta,a_bridge,a_coil,a_turn)
    return numpy.array(sse)


def main():

    data = read_file(sys.argv[1],l=[])
    print data
    

if __name__ == '__main__':
    main()

