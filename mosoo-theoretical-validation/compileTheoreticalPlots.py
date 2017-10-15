#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from __future__ import division
import numpy as np
import paretofront 
import epsilonindicator
import time
import os
from paretofront import paretofront
from paretofront import paretofront_cao
from epsilonindicator import compute_eps
from epsilonindicator import compute_incr_eps
from epsilonindicator import compute_fast_incr_eps
# flags for execution
isPsi = False
isLoss = False
isInd = False
isTbl = True
isCompressTbl = True
NUM_OBJ = 2;
NUM_INSTANCES = 4;
NUM_NUMERICAL = 4
NUM_DIM = 2

HMAX_DIR = "hmax_data"
NUMERICAL_DIR = "numerical_data"
THEORY_DIR = "theoretical_data"
HMAX_DIR = "hmax_data"
PSI_DIR = "psi_data"
PF_DIR = "pf_data"
PF_FILE = "pf-n%d-mosoo-i%d.txt"
EX_PF_FILE = "ext-front-n%d-mosoo-i%d.txt"
PSI_FILE = "psi-n%d-mosoo.txt"
NUMERICAL_FILE = "resultfile-f1-d%d-r%d-fx.txt"
LOSS_FILE = "loss-n%d-mosoo-i%d.txt"
IND_FILE = "ind-n%d-mosoo-i%d.txt"
THRTCL_FILE = "bounds-n%d-mosoo-i%d.txt"
HMAX_FILE = "hmax-n%d-mosoo-i%d.txt"

DIM_RANGE = np.arange(1,3, dtype = int)
ITER_RANGE = np.arange(1,5, dtype = int)

ITER_IDX = 0;
FRST_OBJ_IDX = 1;
SCND_OBJ_IDX = 2;


def compressResultsTable(write = True, isLogScale = False):
    """
    compress the compiled results stored in a table of *.dat extension in the format and overwrites that table
    gnd iter psi loss1 loss2 indicator loss1bound loss2bound indicatorbound
    
    These values are scaled logarthmically are suitable to be used later by pgfplot command
    """
    
    # concatenate the tables after logScaling
    for n in DIM_RANGE:
        psi_results = np.genfromtxt( os.path.join(PSI_DIR, PSI_FILE % (n)),delimiter='\t', comments='%')
        for i in ITER_RANGE:
            print('Compressing tables for dim %d , instance %d' % (n,i))
            print(' 1. Data Table ')
            results = np.genfromtxt( 'table-psi%f-n%d-i%d.dat' % (psi_results[i-1],n, i),delimiter='\t', skip_header=1)
            cmprsd_results  = results.copy()
            # hmax is increasing while other measures are decreasing, hence the different treatment
            row_diff = np.sum(-results[1:-1,1:-1] + results[:-2,1:-1], axis = 1)

            cmprsd_idx = 0
            cmprsd_results[cmprsd_idx,:] = results[0,:]
            for idx, row_val in enumerate(row_diff):
                if row_val > 1e-6:
                    if idx != 0 and prevIdx != idx:
                        cmprsd_idx += 1
                        cmprsd_results[cmprsd_idx,:] = results[idx,:]
                    cmprsd_idx += 1
                    cmprsd_results[cmprsd_idx,:] = results[idx + 1,:]
                    prevIdx = idx + 1
            
            if prevIdx != results.shape[0]-2:
                cmprsd_idx += 1
                cmprsd_results[cmprsd_idx,:] = results[-2,:] 
            cmprsd_idx += 1
            cmprsd_results[cmprsd_idx,:] = results[-1,:]  
            
            data = cmprsd_results[:cmprsd_idx+1,:]
            
            
            np.savetxt('table-psi%f-n%d-i%d.dat' % (psi_results[i-1],n, i), data, fmt='%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f', header='iter psi loss1 loss2 indicator loss1bound loss2bound indicatorbound',comments='')   # use exponential notation
            print(' 2. hmax Table ')
            results = np.genfromtxt( 'table-hmax-n%d-i%d.dat' % (n, i),delimiter='\t', skip_header=1)
            cmprsd_results  = results.copy()
            # hmax is increasing while other measures are decreasing, hence the different treatment
            row_diff = results[1:-1,-1] - results[:-2,-1]

            cmprsd_idx = 0
            cmprsd_results[cmprsd_idx,:] = results[0,:]
            for idx, row_val in enumerate(row_diff):
                if row_val > 1e-6:
                    if idx != 0 and prevIdx != idx:
                        cmprsd_idx += 1
                        cmprsd_results[cmprsd_idx,:] = results[idx,:]
                    cmprsd_idx += 1
                    cmprsd_results[cmprsd_idx,:] = results[idx + 1,:]
                    prevIdx = idx + 1
            
            if prevIdx != results.shape[0]-2:
                cmprsd_idx += 1
                cmprsd_results[cmprsd_idx,:] = results[-2,:] 
            cmprsd_idx += 1
            cmprsd_results[cmprsd_idx,:] = results[-1,:]  
            
            data = cmprsd_results[:cmprsd_idx+1,:]
            
            
            np.savetxt('table-hmax-n%d-i%d.dat' % (n, i), data, fmt='%d\t%d', header='iter hmax',comments='')   # use exponential notation



def generateResultsTable(write = True, isLogScale = False):
    """
    Compile the computed results into a table of *.dat extension in the format
    gnd iter psi loss1 loss2 indicator loss1bound loss2bound indicatorbound
    
    These values are scaled logarthmically are suitable to be used later by pgfplot command
    """
    
    # concatenate the tables after logScaling
    for n in DIM_RANGE:
        for i in ITER_RANGE:
            print('Generating table for dim %d , instance %d' % (n,i))
            loss_results = np.genfromtxt( os.path.join(NUMERICAL_DIR,LOSS_FILE % (n,i)),delimiter='\t', comments='%')
            ind_results = np.genfromtxt( os.path.join(NUMERICAL_DIR,IND_FILE % (n,i)),delimiter='\t', comments='%')
            empirical_results = np.hstack((loss_results, ind_results[:,None]))
            theoretical_results = np.genfromtxt( os.path.join(THEORY_DIR,THRTCL_FILE % (n,i)),delimiter='\t', comments='%')
            psi_results = np.genfromtxt( os.path.join(PSI_DIR,PSI_FILE % (n)),delimiter='\t', comments='%')
            hmax_results = np.genfromtxt( os.path.join(HMAX_DIR,HMAX_FILE % (n,i)),delimiter='\t', comments='%')
            theoretical_results[:,-1] += psi_results[i-1] # add the Psi offset
            
            if isLogScale:
                empirical_results = log(empirical_results);
                theoretical_results = log(theoretical_results);
            # get the data   and clip it to the obtained empirical data (as only non-dominated points are registered) 
            num_points = hmax_results.shape[0]
            # extend empirical_results as only the non-dominated points are recorded:
            empirical_results = np.repeat(empirical_results, [1]*(empirical_results.shape[0]-1) + [1 + num_points- empirical_results.shape[0]], axis=0)
            #gnd = np.zeros((num_points,1))
            iter = np.reshape(np.arange(1,num_points+1),(num_points,1))
            psi = psi_results[i-1] * np.ones((num_points,1))
            #print("Size of iter %d,  psi %d, empirical %d, theorrical %d, hmax %d" %(iter.shape[0],psi.shape[0], empirical_results.shape[0], theoretical_results.shape[0], hmax_results.shape[0]))
            data = np.hstack((iter, psi, empirical_results, theoretical_results))
            hmax_data = np.hstack((iter, hmax_results[:,-1][None,].T))
            
            
            np.savetxt('table-psi%f-n%d-i%d.dat' % (psi_results[i-1],n, i), data, fmt='%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f', header='iter psi loss1 loss2 indicator loss1bound loss2bound indicatorbound',comments='')   # use exponential notation
            np.savetxt('table-hmax-n%d-i%d.dat' % (n, i), hmax_data, fmt='%d\t%d', header='iter hmax',comments='')   # use exponential notation

            
            
def computeLoss(write=True):
    """Compute the loss as a function of the number of iterations for all the files"""
    for n in DIM_RANGE:
        for i in ITER_RANGE:
            print('Processing file: %s' %( NUMERICAL_FILE % (n,i)))
            R = np.genfromtxt( os.path.join(NUMERICAL_DIR,NUMERICAL_FILE % (n,i)),delimiter='\t', comments='%')
            num_iterations = int(np.max(R[:,ITER_IDX]))
            L = np.zeros((num_iterations,NUM_OBJ))
            for j in xrange(num_iterations):
                #print('\t iteration %d' % j)
                idx = np.where(R[:,ITER_IDX]==j+1);
                if not idx[0].size:
                    L[j,0] = L[j-1,0];
                    L[j,1] = L[j-1,1];
                else:
                    lastIterIdx = int(np.max(idx));
                    L[j,0] = np.min(R[:lastIterIdx+1,FRST_OBJ_IDX])
                    L[j,1] = np.min(R[:lastIterIdx+1,SCND_OBJ_IDX])
            
            if (write):
                np.savetxt(os.path.join(NUMERICAL_DIR,LOSS_FILE % (n,i)), L, delimiter='\t')
        
def computeIndicator(write=True):
    """Compute the loss as a function of the number of iterations for all the files"""
    for n in DIM_RANGE:
        for i in ITER_RANGE:
            print('Processing file: %s' %( NUMERICAL_FILE % (n,i)))
            P = np.genfromtxt( os.path.join(PF_DIR,PF_FILE % (n,i)),delimiter='\t')
            EP = np.genfromtxt( os.path.join(PF_DIR,EX_PF_FILE % (n,i)),delimiter='\t')
            # subsample it for ease
            P = P[::5,:]
            P = np.vstack((P,EP))
            R = np.genfromtxt( os.path.join(NUMERICAL_DIR,NUMERICAL_FILE % (n,i)),delimiter='\t', comments='%')
            epsilonVector = compute_fast_incr_eps(R[:,FRST_OBJ_IDX:SCND_OBJ_IDX+1],P)
            num_iterations = int(np.max(R[:,ITER_IDX]))
            I = np.zeros(num_iterations)
            for j in xrange(num_iterations):
                #print("Iteration %d" % j)
                idx = np.where(R[:,ITER_IDX]==j+1);
                if not idx[0].size:
                    I[j] = I[j-1]
                else:
                    lastIterIdx = int(np.max(np.where(R[:,ITER_IDX]==j+1)))
                    #I[j] = computeEps(R[:lastIterIdx+1,FRST_OBJ_IDX:SCND_OBJ_IDX+1],P)
                    I[j] = epsilonVector[lastIterIdx]
            
            if (write):
                np.savetxt(os.path.join(NUMERICAL_DIR,IND_FILE % (n,i)), I, delimiter='\t')       
    


def computePsi(write= True):
    """ Compute the conflict dimension of the problems whose pareto front and
            Pareto front extrema are saved in the files PF_FILE, EX_PF_FILE, respectively.
            NUM_INSTANCES is the number of files (problems) to be processed
    """
    for n in DIM_RANGE:
        Psi = [1] * NUM_INSTANCES;
        for i in ITER_RANGE:
            pf = np.genfromtxt(os.path.join(PF_DIR,PF_FILE % (n, i)), delimiter='\t')
            epf = np.genfromtxt(os.path.join(PF_DIR,EX_PF_FILE % (n, i)), delimiter='\t')
            Psi[i-1]= compute_eps(epf, pf);
            print("Psi for instance %d is %f" % (i, Psi[i-1]));
        if write:
            np.savetxt(os.path.join(PSI_DIR,PSI_FILE % n), Psi, delimiter='\t')
    

# Thanks to Do-Thanh Tran and Dimo for the two functions below:
def computeEps(dataSet, refSet, method="additive"):
    """All objectives are subject to minimization.
    
    dataSet: Is a 2D numpy array (nData x nObj) containing the output of an MO run
             Should be nondominated points.
    
    refSet:  Is a 2D numpy array (nRef x nObj) containing the reference set
             with which the indicator is calculated.
    
    method:  Specifies the indicator type to use: "additive" or "multiplicative".
             Default is "additive" (i.e. the additive epsilon indicator).
    """
    
    if method == "additive":
        value = -np.inf
        additiveEpsilon = True
    elif method == "multiplicative":
        value = 0
        additiveEpsilon = False
    else:
        print "Wrong 'method' parameter!"
        quit()
    
    nData, nObj = dataSet.shape
    nRef, dimSet = refSet.shape
    if nObj < 1 or nObj != dimSet:
        print "Number of objectives mismatched."
        print "Dimensions of dataSet and refSet must be equal and > 1"
        quit()
    
    # For any interpretation of the epsilon indicator,
    # the outer loop must be over the points in refSet.
    for i in xrange(nRef):
        for j in xrange(nData):
            eps_j = epsilonDominance(dataSet[j,:], refSet[i,:], nObj, additiveEpsilon)
            if j == 0:
                eps_i = eps_j   # epsilon from dataSet to point i in refSet
            elif eps_i > eps_j: # takes the SMALLEST epsilon over all points of dataSet to point i in refSet
                eps_i = eps_j
        
        if i == 0:
            value = eps_i   # epsilon from dataSet to refSet
        elif value < eps_i: # takes the LARGEST epsilon over all points of refSet
            value = eps_i   # thereby when transformed by a factor ${value}, the whole refSet is weakly dominated by dataSet
    
    return value


def epsilonDominance(dataPoint, refPoint, nObj, additive=True):
    for k in xrange(nObj):
        if additive: # additive epsilon
            eps_k = dataPoint[k] - refPoint[k] # distance in this k-th objective
        else: # multiplicative
            if (refPoint[k] < 0 and dataPoint[k] > 0) or \
               (refPoint[k] > 0 and dataPoint[k] < 0) or \
               (refPoint[k] == 0 and dataPoint[k] == 0):
                print "dataPoint and refPoint have to be > 0"
                quit()
            eps_k = dataPoint[k] / refPoint[k]
        if k == 0:
            epsilon = eps_k   # epsilon from point j in dataSet to point i in refSet
        elif epsilon < eps_k: # takes the LARGEST over all the objectives
            epsilon = eps_k
    return epsilon
	
	
	

	
if __name__ == "__main__":
    if isPsi:
        print("Computing Psi for the problems specified in me :)")
        computePsi();
        print("I am done!")
    if isLoss:
        print("Computing Loss for the problems specified in me :)")
        computeLoss();
        print("I am done!")
    if isInd:
        print("Computing additive indicator for the problems specified in me :)")
        computeIndicator();
        print("I am done!")
        
    if isTbl:
        print("Generating tables for the problems specified in me :)")
        generateResultsTable()
        print("I am done!")
        
    if isCompressTbl:
        print("Removing redundancy from tables :)")
        compressResultsTable()
        print("I am done!")       