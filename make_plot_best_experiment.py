#! /usr/bin/env python

from matplotlib import pyplot

def loaddata(f):
    fo = open(f,'r')
    name = {}
    for l in fo:
        a = l.split()
        name[a[0].upper()] = float(a[1])
    return name

def loadbest(f):
    fo = open(f,'r')
    name = {}
    for l in fo:
        a = l.split()
        name[a[0].upper()] = float(a[2])
    return name


dictio_exp = loaddata('data_paper.txt')
dictio_sim = loadbest('Results_Best.txt')
norm_exp = dictio_exp['ALA']
norm_sim = dictio_sim['ALA']
#make absolute terms, best results are 10 times smaller than they should
norm_exp = 1
norm_sim = 0.1

aa = {'ALA':'A','VAL':'V','LEU':'L','ILE':'I','PHE':'F','TRP':'W','MET':'M','PRO':'P','ASP':'D','GLU':'E','GLY':'G','SER':'S','THR':'T','CYS':'C','TYR':'Y','ASN':'N','GLN':'Q','LYS':'K','ARG':'R','HIS':'H'}

pyplot.figure()
ax = pyplot.subplot(111)
for a in aa:
    try:
        dictio_exp[a] /= norm_exp
        dictio_sim[a] /= norm_sim
        pyplot.plot(dictio_exp[a],dictio_sim[a],'ko',markersize=8)
        print dictio_exp[a],dictio_sim[a],aa[a]
        #print a,dictio_sim[a]-dictio_sim['ALA']
        ax.text(dictio_exp[a],dictio_sim[a]-0.08,aa[a],size=10,color='black')
    except:
        pass

pyplot.plot(range(3),range(3),'k-')
x = [0,1,2,3,4,5]
y = [0,0.5,1,1.5,2,2.5]
pyplot.plot(x,y,'k--')
pyplot.plot(y,x,'k--')
pyplot.xlim(0,2.0)
pyplot.ylim(0,2.0)
#ax.set_xlim(0,2.0)
#ax.set_ylim(0,2.8)
#ax.set_xlim(0,2.0)
#ax.set_ylim(0,2.0)
ax.set_aspect(1)
pyplot.savefig('experimental.pdf',format='pdf')
