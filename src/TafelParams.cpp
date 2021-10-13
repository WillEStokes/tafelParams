#include <iostream>
#include <array>
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h> // for pow
#include <numeric>
#include <algorithm>
#include <chrono> // for timing
#include <unistd.h> // for sleep
#include <ctime>
#include "PolynomialRegression.h"
#include "SolveQuadratic.h"
#include "SolveCubic.h"
#include "TafelParams.h"

// Function declerations
void readFileWrapper(char* file_name, double data[2817]);
void readFile(std::string file_name, std::vector<double>& data);
std::vector<double> deriv(std::vector<double>& data);
std::vector<double> fastSmooth(std::vector<double>& yData, unsigned int smoothWidth);
std::vector<double> forEachVector(const std::vector<double>& values, const std::vector<double>& modifier, double(*func)(double, double));
std::vector<double> forEachScalar(const std::vector<double>& values, double modifier, double(*func)(double, double));
std::vector<double> polynomial(std::vector<double>& xData, std::vector<double>& coeffs);
std::vector<double> linspace(double startVal, double endVal, int elements);
std::vector<double> logspace(double startVal, double endVal, int elements);
int inRange(double val, double lowerLim, double upperLim);
int getTafelParams(std::vector<double>& xDataAll, 
    std::vector<double>& yDataAll, 
    const double currentThreshold, 
    const double potentialThreshold, 
    const double cathodicCorrCoeff, 
    double*& cathodicConstant, 
    double*& anodicConstant, 
    double*& iCorr, 
    double*& eCorr, 
    std::array<std::array<double, 2>, 2>& cathodicFitCoords, 
    std::array<std::array<double, 2>, 2>& anodicFitCoords,
    std::array<double, 2>& anodicTraceCoeffs,
    std::array<double, 2>& anodicFitCoeffs );
int getTafelParamsWrapper(double* xData, 
    double* yData, 
    int numRows, 
    const double currentThreshold, 
    const double potentialThreshold, 
    const double cathodicCorrCoeff, 
    double* cathodicConstant, 
    double* anodicConstant, 
    double* iCorr, 
    double* eCorr, 
    double cathodicFitCoords[4], 
    double anodicFitCoords[4],
    double anodicTraceCoeffs[2],
    double anodicFitCoeffs[2] );
double anodicFitTest(double* xDataDbL, 
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
double cathodicFitTest(double* xDataDbL, 
    double* yDataDbL, 
    int numRows,  
    double cathodicCorrCoeff,
    double* eCorr,
    int* cathodicElements,
    double xCathodicFit[200],
    double yCathodicFit[200],
    double corrCoeffVals[200],
    double fitConstants[2]);
int leastSquares(std::vector<double> xData,
    std::vector<double> yData,
    int elements,
    int linearFit,
    std::pair<double, double>& constants,
    std::vector<double>& xFit,
    std::vector<double>& yFit);
double correlationCoefficient(std::vector<double> xData,
    std::vector<double> yData,
    int elements);

int main()
{
    std::vector<double> xDataAll;
    std::vector<double> yDataAll;
    readFile("C:/Users/menwst/Documents/CPP/tafelParams/data/xData", xDataAll);
    readFile("C:/Users/menwst/Documents/CPP/tafelParams/data/yData", yDataAll);

    double dummy0 = 0;
    double dummy1 = 0;
    double dummy2 = 0;
    double dummy3 = 0;
    const double currentThreshold = -0.4;
    const double potentialThreshold = 1.1;
    const double cathodicCorrCoeff = 0.35;
    double* cathodicConstant = &dummy0;
    double* anodicConstant = &dummy1;
    double* eCorr = &dummy2;
    double* iCorr = &dummy3;
    std::array<std::array<double, 2>, 2> cathodicFitCoords = {0, 0, 0, 0};
    std::array<std::array<double, 2>, 2> anodicFitCoords = {0, 0, 0, 0};
    std::array<double, 2> anodicTraceCoeffs = {0, 0};
    std::array<double, 2> anodicFitCoeffs = {0, 0};

    int error = getTafelParams(xDataAll, yDataAll, currentThreshold, potentialThreshold, cathodicCorrCoeff, cathodicConstant, anodicConstant, iCorr, eCorr, cathodicFitCoords, anodicFitCoords, anodicTraceCoeffs, anodicFitCoeffs);

    std::cout << "error: " << error << std::endl;
    std::cout << "cathodic constant: " << *cathodicConstant << std::endl;
    std::cout << "anodic constant: " << *anodicConstant << std::endl;
    std::cout << "iCorr: " << *iCorr << std::endl;
    std::cout << "eCorr: " << *eCorr << std::endl;
    std::cout << "anodicFitCoords: " << anodicFitCoords[0][0] << ", " << anodicFitCoords[1][0] << ", " << anodicFitCoords[0][1] << ", " << anodicFitCoords[1][1] << std::endl;
    std::cout << "cathodicFitCoords: " << cathodicFitCoords[0][0] << ", " << cathodicFitCoords[1][0] << ", " << cathodicFitCoords[0][1] << ", " << cathodicFitCoords[1][1] << std::endl;
    std::cin.get();
}

int getTafelParamsWrapper(double* xData, 
    double* yData, 
    int numRows, 
    const double currentThreshold, 
    const double potentialThreshold, 
    const double cathodicCorrCoeff, 
    double* cathodicConstant, 
    double* anodicConstant, 
    double* iCorr, 
    double* eCorr, 
    double cathodicFitCoords[4], 
    double anodicFitCoords[4],
    double anodicTraceCoeffs[2],
    double anodicFitCoeffs[2] )
{
    std::vector<double> xDataVector(xData, xData + numRows);
    std::vector<double> yDataVector(yData, yData + numRows);

    std::array<std::array<double, 2>, 2> cathodicFitCoordsArr = {0, 0, 0, 0};
    std::array<std::array<double, 2>, 2> anodicFitCoordsArr = {0, 0, 0, 0};
    std::array<double, 2> anodicTraceCoeffsArr = {0, 0};
    std::array<double, 2> anodicFitCoeffsArr = {0, 0};
    
    int error = getTafelParams(xDataVector, yDataVector, currentThreshold, potentialThreshold, cathodicCorrCoeff, cathodicConstant, anodicConstant, iCorr, eCorr, cathodicFitCoordsArr, anodicFitCoordsArr, anodicTraceCoeffsArr, anodicFitCoeffsArr);

    int it = 0;
    for (int i = 0; i < 2; i++)
    {
        anodicTraceCoeffs[i] = anodicTraceCoeffsArr[i];
        anodicFitCoeffs[i] = anodicFitCoeffsArr[i];
        for (int j = 0; j < 2; j++)
        {
            cathodicFitCoords[it] = cathodicFitCoordsArr[i][j];
            anodicFitCoords[it] = anodicFitCoordsArr[i][j];
            it = it + 1;
        }
    }

    return error;
}

int getTafelParams(std::vector<double>& xDataAll, 
    std::vector<double>& yDataAll, 
    const double currentThreshold, 
    const double potentialThreshold, 
    const double cathodicCorrCoeff, 
    double*& cathodicConstant, 
    double*& anodicConstant, 
    double*& iCorr, 
    double*& eCorr, 
    std::array<std::array<double, 2>, 2>& cathodicFitCoords, 
    std::array<std::array<double, 2>, 2>& anodicFitCoords,
    std::array<double, 2>& anodicTraceCoeffs,
    std::array<double, 2>& anodicFitCoeffs )
{
    std::vector<double> xData;
    std::vector<double> yData;
    for (int i = 0; i < xDataAll.size(); i++)
    {
        if (xDataAll[i] < potentialThreshold && yDataAll[i] > currentThreshold)
        {
            xData.push_back(xDataAll[i]);
            yData.push_back(yDataAll[i]);
        }
    }

    if (xData.size() < 4) {return 1;}
    
    std::vector<double> xDeriv = fastSmooth(xData, 3);
    xDeriv = deriv(xDeriv);

    xDeriv.erase(xDeriv.begin(), xDeriv.begin() + 2);
    xData.erase(xData.begin(), xData.begin() + 2);
    yData.erase(yData.begin(), yData.begin() + 2);

    for (int i = 0; i < xDeriv.size(); i++)
    {
        if (xDeriv[i] > 0 )
            break;
        if (i + 1 == xDeriv.size())
            return 2;
    }

    double minX = *std::min_element(xData.begin(), xData.end());
    int zeroCrossInd;
    
    for (zeroCrossInd = 0; zeroCrossInd < xData.size(); zeroCrossInd++)
    {
        if (xData[zeroCrossInd] == minX )
        {
            *eCorr = yData[zeroCrossInd];
            break;
        }
    }

    std::vector<double> xCathodic, yCathodic;
    xCathodic.push_back(xData[0]);
    yCathodic.push_back(yData[0]);
    PolynomialRegression<double> polyReg;
    std::vector<double> corrCoeffVect, xCathodicFitVect, yCathodicFitVect;
    std::pair<double, double> constants;
    std::vector<std::pair<double, double>> constantsVect;
    double corrCoeff;
    int j = 2;
    for (int i = 1; i < zeroCrossInd; i++ )
    {
        xCathodic.push_back(xData[i]);
        yCathodic.push_back(yData[i]);
        xCathodicFitVect.clear();
        yCathodicFitVect.clear();

        leastSquares(xCathodic, yCathodic, j, true, constants, xCathodicFitVect, yCathodicFitVect);

        corrCoeff = correlationCoefficient(xCathodicFitVect, yCathodicFitVect, j);

        corrCoeffVect.push_back(corrCoeff);
        constantsVect.push_back(constants);
        j = j + 1;
    }

    int cFitInd;
    for (cFitInd = corrCoeffVect.size() - 1; cFitInd >= 2; cFitInd-- )
    {
        if (cFitInd == 2)
            return 1;
        if (abs(corrCoeffVect[cFitInd]) >= cathodicCorrCoeff)
            break;
    }

    xCathodicFitVect.clear();
    yCathodicFitVect.clear();
    *iCorr = exp((*eCorr - constantsVect[cFitInd].first) / constantsVect[cFitInd].second);
    std::vector<double> xFit = linspace(*iCorr, *std::max_element(xCathodic.begin(), xCathodic.end()), 2);
    std::vector<double> yFit;
    for (int i = 0; i < 2; i++)
        yFit.push_back((constantsVect[cFitInd].first + constantsVect[cFitInd].second * log(xFit[i])));

    cathodicFitCoords = {xFit[0], xFit[1], yFit[0], yFit[1]};
    *cathodicConstant = (yFit[1] - yFit[0]) / (xFit[1] - xFit[0]);

    int anodicElements = xData.size() - zeroCrossInd;

    std::vector<double> xAnodic;
    std::vector<double> yAnodic;
    std::copy(xData.begin() + zeroCrossInd, xData.end(), back_inserter(xAnodic));
    std::copy(yData.begin() + zeroCrossInd, yData.end(), back_inserter(yAnodic));

    std::vector<double> xAnodicTraceVect, yAnodicTraceVect, xAnodicFitVect, yAnodicFitVect;
    std::pair<double, double> traceConstants, fitConstants;
    leastSquares(xAnodic, yAnodic, anodicElements, false, traceConstants, xAnodicTraceVect, yAnodicTraceVect);

    anodicTraceCoeffs[0] = traceConstants.first;
    anodicTraceCoeffs[1] = traceConstants.second;

    double xUpper = *std::max_element(xAnodic.begin(), xAnodic.end());
    std::vector<double> xAnodicGuess = logspace(*iCorr, xUpper, 50);

    double step = 0.00005;
    double y2 = *eCorr - step;  
    double maxX = *std::max_element(xAnodic.begin(), xAnodic.end());
    std::vector<double> yAnodicGuess = linspace(*eCorr, y2, 50);
    int intersects = false;
    while (!intersects)
    {
        y2 = y2 + step;
        if (y2 >= 0) {
            anodicFitCoeffs[0] = fitConstants.first;
            anodicFitCoeffs[1] = fitConstants.second;
            return 4; }

        yAnodicGuess.clear();
        xAnodicFitVect.clear();
        yAnodicFitVect.clear();
        yAnodicGuess = linspace(*eCorr, y2, 50);
        leastSquares(xAnodicGuess, yAnodicGuess, 50, true, fitConstants, xAnodicFitVect, yAnodicFitVect);

        for (int i = 0; i < 50; i++){
            if (traceConstants.first * pow(traceConstants.second, xAnodicFitVect[i]) < yAnodicFitVect[i])
                intersects = true;
        }
    }

    anodicFitCoeffs[0] = fitConstants.first;
    anodicFitCoeffs[1] = fitConstants.second;

    *anodicConstant = (maxX - *eCorr) / (*iCorr - y2);
    anodicFitCoords = {*iCorr, maxX, *eCorr, y2};

    return 0;
}

double cathodicFitTest(double* xDataDbL, 
    double* yDataDbL, 
    int numRows,  
    double cathodicCorrCoeff,
    double* eCorr,
    int* cathodicElements,
    double xCathodicFit[200],
    double yCathodicFit[200],
    double corrCoeffVals[200],
    double fitConstants[2])
{
    std::vector<double> xData(xDataDbL, xDataDbL + numRows);
    std::vector<double> yData(yDataDbL, yDataDbL + numRows);
    
    std::vector<double> xDeriv = fastSmooth(xData, 3);
    xDeriv = deriv(xDeriv);

    xDeriv.erase(xDeriv.begin(), xDeriv.begin() + 2);
    xData.erase(xData.begin(), xData.begin() + 2);
    yData.erase(yData.begin(), yData.begin() + 2);

    for (int i = 0; i < xDeriv.size(); i++)
    {
        if (xDeriv[i] > 0 )
            break;
        if (i + 1 == xDeriv.size())
            return 2;
    }

    double minX = *std::min_element(xData.begin(), xData.end());
    int zeroCrossInd;
    
    for (zeroCrossInd = 0; zeroCrossInd < xData.size(); zeroCrossInd++)
    {
        if (xData[zeroCrossInd] == minX )
        {
            *eCorr = yData[zeroCrossInd];
            break;
        }
    }

    std::vector<double> xCathodic, yCathodic;
    xCathodic.push_back(xData[0]);
    yCathodic.push_back(yData[0]);
    PolynomialRegression<double> polyReg;
    std::vector<double> corrCoeffVect, xCathodicFitVect, yCathodicFitVect;
    std::pair<double, double> constants;
    std::vector<std::pair<double, double>> constantsVect;
    double corrCoeff;
    int j = 2;
    for (int i = 1; i < zeroCrossInd; i++ )
    {
        xCathodic.push_back(xData[i]);
        yCathodic.push_back(yData[i]);
        xCathodicFitVect.clear();
        yCathodicFitVect.clear();

        leastSquares(xCathodic, yCathodic, j, true, constants, xCathodicFitVect, yCathodicFitVect);

        corrCoeff = correlationCoefficient(xCathodicFitVect, yCathodicFitVect, j);

        corrCoeffVect.push_back(corrCoeff);
        constantsVect.push_back(constants);
        j = j + 1;
    }

    for (int i = 0; i < zeroCrossInd - 1; i++)
    {
        corrCoeffVals[i] = corrCoeffVect[i];
    }

    int cFitInd;
    for (cFitInd = corrCoeffVect.size() - 1; cFitInd >= 2; cFitInd-- )
    {
        if (cFitInd == 2)
            return 1;
        if (abs(corrCoeffVect[cFitInd]) >= cathodicCorrCoeff)
            break;
    }

    xCathodicFitVect.clear();
    yCathodicFitVect.clear();
    double xLower = exp((*eCorr - constantsVect[cFitInd].first) / constantsVect[cFitInd].second);
    std::vector<double> xFit = linspace(xLower, *std::max_element(xCathodic.begin(), xCathodic.end()), 2);
    std::vector<double> yFit;
    for (int i = 0; i < 2; i++)
        yFit.push_back((constantsVect[cFitInd].first + constantsVect[cFitInd].second * log(xFit[i])));
    
    fitConstants[0] = constantsVect[cFitInd].first;
    fitConstants[1] = constantsVect[cFitInd].second;
    *cathodicElements = cFitInd;

    for (int i = 0; i < 2; i++)
    {
        xCathodicFit[i] = xFit[i];
        yCathodicFit[i] = yFit[i];
    }

    return 0;
}

double anodicFitTest(double* xDataDbL, 
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
    double anodicFitCoeffs[4])
{
    *intersects = false;

    std::vector<double> xData(xDataDbL, xDataDbL + numRows);
    std::vector<double> yData(yDataDbL, yDataDbL + numRows);
    
    std::vector<double> xDeriv = fastSmooth(xData, 3);
    xDeriv = deriv(xDeriv);

    xDeriv.erase(xDeriv.begin(), xDeriv.begin() + 2);
    xData.erase(xData.begin(), xData.begin() + 2);
    yData.erase(yData.begin(), yData.begin() + 2);

    for (int i = 0; i < xDeriv.size(); i++)
    {
        if (xDeriv[i] > 0 )
            break;
        if (i + 1 == xDeriv.size())
            return 2;
    }

    double minX = *std::min_element(xData.begin(), xData.end());
    int zeroCrossInd;
    
    for (zeroCrossInd = 0; zeroCrossInd < xData.size(); zeroCrossInd++)
    {
        if (xData[zeroCrossInd] == minX )
        {
            *eCorr = yData[zeroCrossInd];
            break;
        }
    }

    *anodicElements = xData.size() - zeroCrossInd;

    std::vector<double> xAnodic;
    std::vector<double> yAnodic;
    std::copy(xData.begin() + zeroCrossInd, xData.end(), back_inserter(xAnodic));
    std::copy(yData.begin() + zeroCrossInd, yData.end(), back_inserter(yAnodic));

    std::vector<double> xAnodicTraceVect, yAnodicTraceVect, xAnodicFitVect, yAnodicFitVect;
    std::pair<double, double> traceConstants, fitConstants;
    leastSquares(xAnodic, yAnodic, *anodicElements, false, traceConstants, xAnodicTraceVect, yAnodicTraceVect);

    double xUpper = *std::max_element(xAnodic.begin(), xAnodic.end());
    std::vector<double> xAnodicGuess = logspace(iCorr, xUpper, 50);
    std::vector<double> yAnodicGuess = linspace(*eCorr, y2, 50);
    leastSquares(xAnodicGuess, yAnodicGuess, 50, true, fitConstants, xAnodicFitVect, yAnodicFitVect);

    for (int i = 0; i < 50; i++)
    {
        if (traceConstants.first * pow(traceConstants.second, xAnodicFitVect[i]) < yAnodicFitVect[i]) {
            *intersects = true;
            break; }
    }

    for (int i = 0; i < *anodicElements; i++)
    {
        xAnodicTrace[i] = xAnodicTraceVect[i];
        yAnodicTrace[i] = yAnodicTraceVect[i];
    }

    for (int i = 0; i < 50; i++)
    {
        xAnodicFit[i] = xAnodicFitVect[i];
        yAnodicFit[i] = yAnodicFitVect[i];
    }

    anodicTraceCoeffs[0] = traceConstants.first;
    anodicTraceCoeffs[1] = traceConstants.second;

    anodicFitCoeffs[0] = fitConstants.first;
    anodicFitCoeffs[1] = fitConstants.second;

    return 0;
}

int leastSquares(std::vector<double> xData,
    std::vector<double> yData,
    int elements,
    int linearFit,
    std::pair<double, double>& constants,
    std::vector<double>& xFit,
    std::vector<double>& yFit)
{
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, sumLogX = 0, sumLogX2 = 0, sumLogXY = 0, sumLogY = 0, sumXLogY = 0, a, b, A, B;

    std::cout << "------LEASTSQUARES OUTPUT-------" << std::endl;

    for (int i = 0; i < elements; i++)
    {
        sumX = sumX + xData[i];
        sumX2 = sumX2 + xData[i] * xData[i];
        sumXY = sumXY + xData[i] * yData[i];
        sumY = sumY + yData[i];
        sumLogX = sumLogX + log(xData[i]);
        sumLogX2 = sumLogX2 + log(xData[i]) * log(xData[i]);
        sumLogXY = sumLogXY + log(xData[i]) * yData[i];
        sumLogY = sumLogY + log(abs(yData[i]));
        sumXLogY = sumXLogY + xData[i] * log(abs(yData[i]));
        std::cout << "xCathodic: " << xData[i] << "; yCathodic: " << yData[i] << std::endl;
    }

    std::cout << std::endl;
    std::cout << "sumX: " << sumX << "; sumX2: " << sumX2 << "; sumY: " << sumY << "; sumXY: " << sumXY << std::endl;
    std::cout << std::endl;

    if (linearFit)
    {
        B = (elements * sumLogXY - sumLogX * sumY) / (elements * sumLogX2 - sumLogX * sumLogX);
        A = (sumY - B * sumLogX) / elements;
        constants.first = A;
        constants.second = B;
    }
    else
    {
        B = (elements * sumXLogY - sumX * sumLogY) / (elements * sumX2 - sumX * sumX);
        A = (sumLogY - B * sumX) / elements;
        a = -exp(A);
        b = exp(B);
        constants.first = a;
        constants.second = b;
    }

    std::cout << "A: " << A << "; B: " << B << "; a: " << a << "; b: " << b << std::endl;
    std::cout << std::endl;

    xFit = logspace(*std::min_element(xData.begin(), xData.begin() + elements - 1), *std::max_element(xData.begin(), xData.begin() + elements - 1), elements);
    for (int i = 0; i < elements; i++)
    {
        if (linearFit) {
            yFit.push_back((A + B * log(xFit[i]))); }
        else {
            yFit.push_back(a * pow(b, xFit[i])); }
        std::cout << "xFit: " << xFit[i] << "; yFit: " << yFit[i] << std::endl;
    }

    std::cout << std::endl;

    return 0;
}

double correlationCoefficient(std::vector<double> xData,
    std::vector<double> yData,
    int elements)
{
    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double squareSum_X = 0, squareSum_Y = 0;

    std::cout << "------CORRELATIONCOEFFICIENT OUTPUT-------" << std::endl;
  
    for (int i = 0; i < elements; i++)
    {
        sum_X = sum_X + xData[i];
        sum_Y = sum_Y + abs(yData[i]);
        sum_XY = sum_XY + xData[i] * abs(yData[i]);
        squareSum_X = squareSum_X + xData[i] * xData[i];
        squareSum_Y = squareSum_Y + abs(yData[i]) * abs(yData[i]);
        std::cout << "xFit: " << xData[i] << "; yFit: " << yData[i] << std::endl;
    }

    std::cout << std::endl;
    std::cout << "sum_X: " << sum_X << "; sum_Y: " << sum_Y << "; sum_XY: " << sum_XY << "; squareSum_X: " << squareSum_X << "; squareSum_Y: " << squareSum_Y << std::endl;
    std::cout << std::endl;

    double corrCoeff = (double)(elements * sum_XY - sum_X * sum_Y) / sqrt((elements * squareSum_X - sum_X * sum_X) * (elements * squareSum_Y - sum_Y * sum_Y));

    std::cout << "corrCoeff: " << corrCoeff << std::endl;
    std::cout << std::endl;
  
    return corrCoeff;
}

int inRange(double val, double lowerLim, double upperLim)
{
    if (val > lowerLim && val < upperLim)
        return true;
    else
        return false;
}

std::vector<double> linspace(double startVal, double endVal, int elements)
{
    std::vector<double> linspaceVector;

    double delta = (endVal - startVal) / (elements - 1);

    for(int i = 0; i < elements - 1; ++i)
    {
      linspaceVector.push_back(startVal + delta * i);
    }
    linspaceVector.push_back(endVal);

    return linspaceVector;
}

std::vector<double> logspace(double startVal, double endVal, int elements)
{
    std::vector<double> logspaceVector;

    double delta = (log10(endVal) - log10(startVal)) / (elements - 1);

    for (int i = 0; i < elements; i++)
    {
        logspaceVector.push_back(pow(10, log10(startVal) + delta * i));
    }

    return logspaceVector;
}

std::vector<double> polynomial(std::vector<double>& xData, std::vector<double>& coeffs)
{
    std::vector<double> yFit;
    for (int i = 0; i < xData.size(); i++)
    {
        if (coeffs.size() == 2) // 1st order poly
            yFit.push_back(coeffs[0] + coeffs[1] * xData[i]);
        else if (coeffs.size() == 3) // 2nd order poly
            yFit.push_back(coeffs[0] + coeffs[1] * xData[i] + coeffs[2] * xData[i] * xData[i]);
        else if (coeffs.size() == 4) // 3rd order poly
            yFit.push_back(coeffs[0] + coeffs[1] * xData[i] + coeffs[2] * xData[i] * xData[i] + coeffs[3] * xData[i] * xData[i] * xData[i]);
    }

    return yFit;
}

void readFile(std::string file_name, std::vector<double>& data)
{
    std::ifstream file;
    file.open(file_name);

    std::string column_one;

    while(getline(file, column_one))
    {
        data.push_back(atof(column_one.c_str()));
    }
}

std::vector<double> deriv(std::vector<double>& data)
{
    int numElements = data.size();
    std::vector<double> deriv;
    deriv.push_back(data[1] - data[0]);
    for(int j = 1; j < numElements - 1; j++)
    {
        deriv.push_back((data[j + 1] - data[j - 1]) / 2);
    }
    deriv.push_back(data[numElements - 1] - data[numElements - 2]);

    return deriv;
}

std::vector<double> fastSmooth(std::vector<double>& yData, unsigned int smoothWidth)
{
    std::vector<double> smooth_yData;
    std::vector<double> temp;
    double sumPoints = std::accumulate(yData.begin(), yData.begin() + smoothWidth, 0.0);
    int halfWidth = round(smoothWidth / 2);
    int numPoints = yData.size();
    for (int i = 0; i < halfWidth; i++)
    {
        temp.push_back(0);
    }
    for (int i = 0; i < numPoints - smoothWidth; i++)
    {
        temp.push_back(sumPoints);
        sumPoints = sumPoints - yData[i];
        sumPoints = sumPoints + yData[i + smoothWidth];
    }

    for (int i = 0; i < numPoints; i++)
    {
        smooth_yData.push_back(temp[i] / smoothWidth);
    }

    return smooth_yData;
}

std::vector<double> forEachVector(const std::vector<double>& values, const std::vector<double>& modifier, double(*func)(double, double))
{
    std::vector<double> modifiedArray;
    for (int i = 0; i < values.size(); i++)
        modifiedArray.push_back(func(values[i], modifier[i]));

    return modifiedArray;
}

std::vector<double> forEachScalar(const std::vector<double>& values, double modifier, double(*func)(double, double))
{
    std::vector<double> modifiedArray;
    for (int i = 0; i < values.size(); i++)
        modifiedArray.push_back(func(values[i], modifier));

    return modifiedArray;
}

