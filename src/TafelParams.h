#ifndef TAFELPARAMS_H
#define TAFELPARAMS_H

#define ELEMENTS 1000

/* Include Files */

/* Function Declarations */
extern "C" int getTafelParamsWrapper(double* xData, 
    double* yData, 
    int numRows,
    const double potentialThreshold, 
    const double cathodicRSquared, 
    double cathodicConstant, 
    double anodicConstant, 
    double iCorr, 
    double eCorr, 
    double cathodicFitCoords[4], 
    double anodicFitCoords[4] );

#endif
