#ifndef TAFELPARAMS_H
#define TAFELPARAMS_H

/* Include Files */

/* Function Declarations */
extern "C" int getTafelParamsWrapper(double* xData, 
    double* yData, 
    int numRows, 
    const double currentThreshold, 
    const double potentialThreshold, 
    const double cathodicCorrCoeff, 
    const int fitType, 
    double* cathodicConstant, 
    double* anodicConstant, 
    double* iCorr, 
    double* eCorr, 
    double cathodicFitCoords[4], 
    double anodicFitCoords[4],
    double anodicTraceCoeffs[4],
    double anodicFitCoeffs[2] );;
extern "C" double anodicFitTest(double* xDataDbL, 
    double* yDataDbL, 
    int numRows,  
    double y2,
    double iCorr,
    double* eCorr,
    int* intersects,
    int* anodicElements,
    double xAnodicTrace[250],
    double yAnodicTrace[250],
    double xAnodicFit[50],
    double yAnodicFit[50],
    double solutions[6],
    double anodicTraceCoeffs[4],
    double anodicFitCoeffs[4] );
extern "C" double cathodicFitTest(double* xDataDbL, 
    double* yDataDbL, 
    int numRows,  
    double cathodicCorrCoeff,
    double* eCorr,
    int* cathodicElements,
    double xCathodicFit[200],
    double yCathodicFit[200],
    double corrCoeffVals[200],
    double fitConstants[2]);

#endif
