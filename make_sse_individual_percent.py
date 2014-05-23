#! /usr/bin/env python

import numpy
import sys

prot = sys.argv[2] 
ff = sys.argv[3]

data = numpy.genfromtxt(sys.argv[1],dtype='str').view(numpy.chararray)
l_chain = len(data[0])
l_sim = len(data)

alpha = []
beta = []
for i in range(l_chain):
    tot_hel = numpy.float(numpy.sum(data[:,i].count('H'))) + numpy.float(numpy.sum(data[:,i].count('G'))) 
    alpha.append(tot_hel/l_sim*100)
    beta.append(numpy.float(numpy.sum(data[:,i].count('E')))/l_sim*100)

alpha = numpy.array(alpha)
beta = numpy.array(beta)

numpy.set_printoptions(precision=2)
print prot,ff,"alpha",
for a in alpha:
    print "{} ".format(a),
print
print prot,ff,"beta",
for b in beta:
    print "{} ".format(b),
print
