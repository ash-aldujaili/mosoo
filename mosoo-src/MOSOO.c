#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
//#include "global.h"
#include "paretofiltering.h"
#include "myarray.h"
#include "wrapbbob.h"
#include "myrandom.h"

// TO DO : complete its implementation
size_t h_max(size_t nVars, size_t nEvals,  size_t totalEvals, size_t minLevel, size_t K, size_t iter) {
  
  
  
  double diffN = log(2*totalEvals- nEvals + 1) / log(K - 1);
  double pnVars = sqrt(nVars * nVars * nVars);
  
  //double logN3;
  //double logN2;
  
  //logN3 =   nVars * 10 * sqrt(diffN * diffN * diffN);
  //logN2 =   log(nEvals + 1) * log(nEvals + 1);
  
  double depth =  diffN;
  return (size_t) (minLevel + depth + pnVars);

}

// MOSOO algorithm
void MOSOO(size_t nVars, size_t nObjs, size_t maxFuncEvals, size_t *func_id, size_t *inst_id, double loBound, double upBound, double pDepth, size_t partitionFactor) {
			   
	//printf("Func id %zu\n",*func_id);
	// MOSOO parameters
	size_t totalNumPoints = maxFuncEvals;
	size_t K = partitionFactor;
	size_t midK = K/2;
	//size_t uniformTreeOffset = (size_t) (log(totalNumPoints) / log(K));// for tree depth;
	
	char isVerbose = 0;
	char isOdd = (char) K % 2 != 0;
	char isMaxDepth = 0;
	size_t  curDim;
	size_t  curLevel;
	size_t  startIdx = 0; // for efficient loop over the variables
	double curScale = 0.;
	double offset = 0.;
	double deltaBound = upBound - loBound;
	// MOSOO data-structures (pointers of pointers)
	struct array1Dp pPoints;
	struct array1Dp pObjs;
	double *objVectors;
	double *pointVectors;
	int *pointLevels;
	// a set of flags for the sample points
	// TO DO : for even K we can replace its item
	bool *isExpandable ; // is it expandable currently
	bool *isV; // is it part of the 
  // timing stuff
  clock_t begin, end;
  double time_spent;
	
	// allocate
	objVectors = malloc(sizeof(double) * nObjs * totalNumPoints);
	pointVectors = malloc(sizeof(double) * nVars * totalNumPoints);
	isExpandable = malloc(sizeof(bool) * totalNumPoints);
	isV = malloc(sizeof(bool) * totalNumPoints);
	pointLevels  = malloc(sizeof(int) * totalNumPoints);	
	array1Dp_construct(&pPoints, totalNumPoints);
  array1Dp_construct(&pObjs, totalNumPoints);
	
	// initialization 
	size_t numPoints = 1;
	size_t pointIdx = 0;
	// to keep track of the present levels (depths) in the tree
	size_t minLevel = 1;
	size_t nExpMinLvl = 0;
	size_t maxLevel = 0;
	size_t iter = 0;
	// point-wise initialization
	for(size_t p=0; p < totalNumPoints; p++) {
	  *array1Dp_element(&pPoints, p) = &(pointVectors[p * nVars]);
	  *array1Dp_element(&pObjs, p) = &(objVectors[p * nObjs]);
	  isExpandable[p] = false;
	  isV[p] = false;
	}
  isExpandable[0] = true;
  isV[0] = true;
	pointLevels[0] = 1;
  // sample the center
	for(size_t j =0; j<nVars; j++)  pointVectors[pointIdx * nVars + j] = 0.5 * deltaBound + loBound;
	mobbob_eval_testLogging_mosoo(pObjs.elements, pPoints.elements, func_id, inst_id, nVars, numPoints, nObjs, pointIdx);
	
	curLevel=1;	   
	// Core of the algorithm
	while (true) {
	  iter = iter + 1;
    // 1. Expand Q nodes sequentially and sample their representative sets
    // a) Expand
    pointIdx = numPoints;
    //curLevel = pointLevels[p];
    curDim = (curLevel - 1) % nVars;
    curScale = deltaBound / pow(K, (curLevel - 1) / nVars + 1);
    
    if (isVerbose) printf("=== A NEW ROUND===\n");
    for(size_t p=startIdx; p < numPoints; p++) if (isExpandable[p]) {
		  // online update minLevel
		  if (pointLevels[p] == minLevel) {
		    nExpMinLvl = nExpMinLvl + 1;
		    if (nExpMinLvl == (size_t) pow(K, (minLevel-1))) minLevel = minLevel + 1;
		  }
      // if K is even expnaded nodes (i.e, samples are never visited again, this is simulated
      // by assigning them depth of 0 while working with the actual tree from depth 1
      if (!isOdd) pointLevels[p] = 0; 
      isExpandable[p] = false;
		  // generate the new nodes
      for (size_t k=0; k < K; k++) {
        if (isOdd && k == midK) { // if it is odd, skip the center point as it will be the same as kids
          pointLevels[p] = curLevel + 1;
          continue;
        }
        // update the dimesnion-wise coordinates and other info
        if (k < midK) offset = - curScale * (midK - k - 0.5 * !isOdd) ;
        else offset = curScale * (k - midK + 0.5 * !isOdd);
		    for(size_t j =0; j<nVars; j++)  {
		      if (j == curDim) pointVectors[pointIdx * nVars + j] = pointVectors[p * nVars + j] + offset;
		      else pointVectors[pointIdx * nVars + j] = pointVectors[p * nVars + j];
		    }
		    // update point-wise info
		    pointLevels[pointIdx] = curLevel + 1;
		    isExpandable[pointIdx] = false;
		    // online update maxLevel
		    if (maxLevel < pointLevels[pointIdx]) maxLevel = pointLevels[pointIdx];
		    // move to the next point
		    pointIdx = pointIdx + 1;
		    // force exit if you run out of budget (but sample before exit)
		    if (pointIdx == totalNumPoints) goto eval_children;
		  }  
		  
		  // Display some info
      if (isVerbose){
        printf("Exp node. %zu: (", p);
        for (size_t j=0; j < nVars- 1; j ++) printf("%f,", pointVectors[p * nVars +j]);
        printf("%f), depth=%zu, splitDim=%zu, Obj. vector:(", pointVectors[p * nVars + nVars  - 1], curLevel, curDim);
        for (size_t j=0; j < nObjs- 1; j ++) printf("%f,", objVectors[p * nObjs +j]);
        printf("%f).\n", objVectors[p * nObjs + nObjs - 1]);
        
      }
    
    }
	  // b) Sample
	  //printf("quick test: start:%zu, end:%zu\n", pointIdx, numPoints + pointIdx);
	  eval_children:
	  mobbob_eval_testLogging_mosoo(pObjs.elements, pPoints.elements, func_id, inst_id, nVars, pointIdx - numPoints, nObjs, numPoints);	
    if (pointIdx == totalNumPoints) break;
    numPoints = pointIdx;
    
    
    
    // 2. Identify potentially optimal nodes at the current level
    // a. Advance to the next level and update V
    curLevel = curLevel + 1;
    // check if the current depth within limits
    isMaxDepth = (curLevel > h_max(nVars, numPoints, totalNumPoints,minLevel, K, iter)) || (curLevel > maxLevel);
    if (isMaxDepth) curLevel = minLevel;
    // speed up the coming loops
    for (size_t p = 0; p < numPoints; p++) if (pointLevels[p] == curLevel) {
      startIdx = p;
      break;
    }

    //size_t numCompared= 0;
    if (isMaxDepth) {
      // set past isV to false and nodes at the current level as true
      for (size_t p = 0; p < numPoints; p++) 
        if (pointLevels[p] == curLevel) {
          isV[p] = true; 
          //numCompared = numCompared + 1;
        }
        else 
          isV[p] = false; 
    }
    else {
      // update isV
      for (size_t p = startIdx; p < numPoints; p++) 
        if (pointLevels[p] == curLevel) {
          isV[p] = true;
          //numCompared = numCompared + 1;
        }
    }
    //printf("Number of nodes compared %zu\n", numCompared);
    
    // b. Filter the nodes in V
    //begin = clock();
    selectiveparetofront(isV, objVectors, numPoints, nObjs);
    //end = clock();
    //time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    //printf("Time spent %f\n", time_spent);
    
    for(size_t p = startIdx; p < numPoints; p++) if (pointLevels[p] == curLevel) isExpandable[p] = isV[p];
    
    // monitor the sampling
    if (isVerbose) sleep(3);
				
	}
			 
	// free memory:
	free(isExpandable);
	free(pointLevels);
	free(objVectors);
	free(pointVectors);
	array1Dp_destruct(&pPoints);
  array1Dp_destruct(&pObjs);  		   
			   
}
