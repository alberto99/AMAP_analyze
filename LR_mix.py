#! /usr/bin/env python

import sys
import numpy
from pyevolve import G1DList, GSimpleGA, Initializators, Mutators, Consts, Selectors
from pyevolve import Scaling



'''
Use the Lifson - Roig Model to find helical propensities of different amino acids.
J. Phys Chem B,2009,113,9004-9015

Define an alpha region: -100 <= phi <= -30
-67 < psi < -7
If a residue belongs to that area is considered helical. Additionally, if residues i-1 and i+1 belong to a helical area, then residue i is considered to belong to a helical segment

The different states are: coil, start/end of helix and within helix and their weights are 1, vi and wi

whithin this model, the partition function is defined as:

    z=(0 0 1)Prod_i=1toN Mi (0 1 1)T

                ( wi vi 0)
    where M =   ( 0  0  1)
                ( vi vi 1)

Given this, we optimize ln L = sum_i Nw,i ln(wi) + sum_i Nv,i ln(vi) - N ln(Z)
in a MC criteria, accept if ln(L) increases and accept with prob Ltrial/L if ln(L) decreases

'''

n_steps = 100
v_ALA = 0.54
w_ALA = 1.36

#After run 10
v_ALA = 0.2528024035614776
w_ALA = 1.3932075952692177
def tovw(mat):
    '''
    make a logical matrix:  
        0 --> coil
        1 --> helical segment
        -1--> helical
    '''

    my_s = mat.shape
    h_mat = numpy.zeros( my_s )

    for i in range(len(mat)):
        if mat[i,0] == 1:
            h_mat[i,0] = -1
        if mat[i,len(mat[i])-1] == 1:
            h_mat[i,len(mat[i])-1] = -1
        for j in range(1,len(mat[i])-1):
            if mat[i,j] == 0:
                continue
            if (mat[i,j-1] == 0) or (mat[i,j+1] == 0):
                h_mat[i,j] = -1
                continue
            h_mat[i,j] = 1
    return h_mat

def make_partition(res,v=10,w=100):

    r1 = [w,v,0]
    r2 = [0,0,1]
    r3 = [v,v,1]
    m  = [r1,r2,r3]
    m = numpy.matrix(m)
    c = numpy.matrix(numpy.identity( 3 ))
    for i in range(res):
        c = c * m

    a =  numpy.matrix([0,1,1])
    b = numpy.matrix([0,0,1])
    c = b* c * a.T

    return numpy.matrix.item(c)

def make_Best_partition(res,v=10,w=100):
    '''Robert Best uses AAXAAAAXAAAAXAA = Ac-(AAXAA)3-NH2
    to determine propensities. First calculate for ALA and based on those
    numbers, calculate for the rest'''
    r1 = [w,v,0]
    r2 = [0,0,1]
    r3 = [v,v,1]
    m  = [r1,r2,r3]
    m = numpy.matrix(m)

    r1_ALA = [w_ALA,v_ALA,0]
    r2_ALA = [0,0,1]
    r3_ALA = [v_ALA,v_ALA,1]
    m_ALA  = [r1_ALA,r2_ALA,r3_ALA]
    m_ALA = numpy.matrix(m_ALA)

    c = numpy.matrix(numpy.identity( 3 ))
    for Aa in 'AAXAAAAXAAAAXAA':
        if "A" in Aa:
            c = c * m_ALA
        if "X" in Aa:
            c = c * m

    a =  numpy.matrix([0,1,1])
    b = numpy.matrix([0,0,1])
    c = b* c * a.T

    return numpy.matrix.item(c)

def scoring_function_factory(Nconf,Nw,Nv,Nw_x,Nv_x,N_res):
    def eval(chromosome):
        v = chromosome[0]
        w = chromosome[1]*10
        Z = make_Best_partition(N_res,v=v,w=w)
        
        lnL = Nw * numpy.log(w_ALA) + Nv * numpy.log(v_ALA) + Nw_x * numpy.log(w) + Nv_x * numpy.log(v) - Nconf * numpy.log(Z)
        '''Will be a negative value, want this number to increase,multiply by -1 and minimize instead of maximize (the algorithm cannot handle negative scores)'''
        return -lnL
    return eval


def main(f_phi,f_psi):
    phi = numpy.loadtxt(f_phi)
    psi = numpy.loadtxt(f_psi)

    l_phi = (phi <= -30.) * (phi > -100.) * 1
    l_psi = (psi <= -7.) * (psi >= -67.) *1

    h_matrix = l_phi*l_psi

    logical_mat = tovw(h_matrix)

    number_conf = len(logical_mat)
    number_helical_w = ((logical_mat > 0) * 1).sum()
    number_helical_v = ((logical_mat < 0) * 1).sum()
    number_helical_w_resx = ((logical_mat[:,2] > 0) * 1).sum() +  ((logical_mat[:,7] > 0) * 1).sum() +  ((logical_mat[:,12] > 0) * 1).sum()
    number_helical_v_resx = ((logical_mat[:,2] < 0) * 1).sum() +  ((logical_mat[:,7] < 0) * 1).sum() +  ((logical_mat[:,12] < 0) * 1).sum()

    number_helical_w -= number_helical_w_resx 
    number_helical_v -= number_helical_v_resx 

    n_res = len(logical_mat[0])
    partition_f = make_Best_partition(len(logical_mat[0]),v=1,w=1)

    print logical_mat
    print number_conf,number_helical_w,number_helical_v,partition_f,number_helical_w_resx,number_helical_v_resx

    eval_func = scoring_function_factory(number_conf,number_helical_w,number_helical_v,number_helical_w_resx,number_helical_v_resx,n_res)

    genome = G1DList.G1DList(2)
    genome.evaluator.set(eval_func)
    genome.setParams(rangemin=0.000001, rangemax=1.0, gauss_sigma=0.01)
    genome.initializator.set(Initializators.G1DListInitializatorReal)
    genome.mutator.set(Mutators.G1DListMutatorRealGaussian)
    ga = GSimpleGA.GSimpleGA(genome)
    ga.setMinimax( Consts.minimaxType["minimize"] )
    ga.setElitism(True)
    #ga.selector.set(Selectors.GRouletteWheel)
    ga.selector.set(Selectors.GRankSelector)
    ga.nGenerations = n_steps
    ga.setMutationRate(0.20)
    ga.evolve(freq_stats=10)
    print ga.bestIndividual()

def eval_func(chromosome):
   score = 0.0
   # iterate over the chromosome
   for value in chromosome:
      if value==0:
         score += 1
   return score





if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2])

