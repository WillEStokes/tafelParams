#include <iostream>
#include <array>
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h> // pow
#include <numeric>
#include <algorithm>
#include <chrono> // for timing
#include <unistd.h> // for sleep
#include <ctime>
#include "PolynomialRegression.h"
#include "SolveQuadratic.h"
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
int inRange(double val, double lowerLim, double upperLim);
int getTafelParams(std::vector<double>& xDataAll, 
    std::vector<double>& yDataAll, 
    const double potentialThreshold, 
    const double cathodicRSquared, 
    double& cathodicConstant, 
    double& anodicConstant, 
    double& iCorr, 
    double& eCorr, 
    std::array<std::array<double, 2>, 2>& cathodicFitCoords, 
    std::array<std::array<double, 2>, 2>& anodicFitCoords );
int getTafelParamsWrapper(double* xData, 
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

int main()
{
    std::vector<double> xDataAll;
    std::vector<double> yDataAll;
    readFile("C:/Users/menwst/Documents/CPP/tafelParams/data/xData", xDataAll);
    readFile("C:/Users/menwst/Documents/CPP/tafelParams/data/yData", yDataAll);

    const double potentialThreshold = 0.01;
    const double cathodicRSquared = 0.9;
    double cathodicConstant = 0;
    double anodicConstant = 0;
    double iCorr = 0;
    double eCorr = 0;
    std::array<std::array<double, 2>, 2> cathodicFitCoords = {0, 0, 0, 0};
    std::array<std::array<double, 2>, 2> anodicFitCoords = {0, 0, 0, 0};

    int error = getTafelParams(xDataAll, yDataAll, potentialThreshold, cathodicRSquared, cathodicConstant, anodicConstant, iCorr, eCorr, cathodicFitCoords, anodicFitCoords);
    
    std::cout << "error: " << error << std::endl;
    std::cout << "cathodic constant: " << cathodicConstant << std::endl;
    std::cout << "anodic constant: " << anodicConstant << std::endl;
    std::cout << "iCorr: " << iCorr << std::endl;
    std::cout << "eCorr: " << eCorr << std::endl;
    std::cout << "anodicFitCoords: " << anodicFitCoords[0][0] << ", " << anodicFitCoords[1][0] << ", " << anodicFitCoords[0][1] << ", " << anodicFitCoords[1][1] << std::endl;
    std::cout << "cathodicFitCoords: " << cathodicFitCoords[0][0] << ", " << cathodicFitCoords[1][0] << ", " << cathodicFitCoords[0][1] << ", " << cathodicFitCoords[1][1] << std::endl;
    std::cin.get();
}

int getTafelParamsWrapper(double* xData, 
    double* yData, 
    int numRows,
    const double potentialThreshold, 
    const double cathodicRSquared, 
    double cathodicConstant, 
    double anodicConstant, 
    double iCorr, 
    double eCorr, 
    double cathodicFitCoords[4], 
    double anodicFitCoords[4] )
{
    std::vector<double> xDataVector;
    std::vector<double> yDataVector;
    for (int i = 0; i < numRows; i++)
    {
        xDataVector.push_back(xData[i]);
        yDataVector.push_back(yData[i]);
    }

    std::array<std::array<double, 2>, 2> cathodicFitCoordsArr = {0, 0, 0, 0};
    std::array<std::array<double, 2>, 2> anodicFitCoordsArr = {0, 0, 0, 0};
    int error = getTafelParams(xDataVector, yDataVector, potentialThreshold, cathodicRSquared, cathodicConstant, anodicConstant, iCorr, eCorr, cathodicFitCoordsArr, anodicFitCoordsArr);

    int it = 0;
    for (int i = 0; i < 2; i++)
    {
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
    const double potentialThreshold, 
    const double cathodicRSquared, 
    double& cathodicConstant, 
    double& anodicConstant, 
    double& iCorr, 
    double& eCorr, 
    std::array<std::array<double, 2>, 2>& cathodicFitCoords, 
    std::array<std::array<double, 2>, 2>& anodicFitCoords )
{
    std::vector<double> xData;
    std::vector<double> yData;
    for (int i = 0; i < xDataAll.size(); i++)
    {
        if (xDataAll[i] < potentialThreshold)
        {
            xData.push_back(xDataAll[i]);
            yData.push_back(yDataAll[i]);
        }
    }

    std::cout << "raw data: " << std::endl;
    for(int i = 0; i < xData.size(); i++ )
    {
        std::cout << xData[i] << ", " << yData[i] << std::endl;
    }
    
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
            return 1;
    }

    double minX = *std::min_element(xData.begin(), xData.end());
    int zeroCrossInd;
    for (zeroCrossInd = 0; zeroCrossInd < xData.size(); zeroCrossInd++)
    {
        if (xData[zeroCrossInd] == minX )
        {
            eCorr = yData[zeroCrossInd];
            break;
        }
    }
    
    std::vector<double> xCathodic;
    std::vector<double> yCathodic;
    std::vector<double> xAnodic;
    std::vector<double> yAnodic;
    std::copy(xData.begin() + zeroCrossInd, xData.end(), back_inserter(xAnodic));
    std::copy(yData.begin() + zeroCrossInd, yData.end(), back_inserter(yAnodic));

    xCathodic.push_back(xData[0]);
    yCathodic.push_back(yData[0]);
    std::vector<double> logXCathodic;
    logXCathodic.push_back(log(xData[0]));
    PolynomialRegression<double> polyReg;
    std::vector<double> cathodicCoeffs;
    std::vector<double> rSquaredArray;
    std::vector<double> ss_res_array;
    std::vector<double> ss_tot_array;
    std::vector<double> yCathodicFit;
    double yMean, ss_res, ss_tot;
    for (int i = 1; i < zeroCrossInd; i++)
    {
        cathodicCoeffs.clear();
        yCathodicFit.clear();
        
        xCathodic.push_back(xData[i]);
        logXCathodic.push_back(log(xData[i]));
        yCathodic.push_back(yData[i]);

        polyReg.fitIt(logXCathodic, yCathodic, 1, cathodicCoeffs );

        yCathodicFit = polynomial(logXCathodic, cathodicCoeffs);

        yMean = std::accumulate(yCathodic.begin(), yCathodic.end(), 0.0) / yCathodic.size();

        ss_res_array = forEachVector(yCathodic, yCathodicFit, [](double a, double b) ->double {return pow(a - b, 2); });
        ss_tot_array = forEachScalar(yCathodic, yMean, [](double a, double b) ->double {return pow(a - b, 2); });
        ss_res = std::accumulate(ss_res_array.begin(), ss_res_array.end(), 0.0);
        ss_tot = std::accumulate(ss_tot_array.begin(), ss_tot_array.end(), 0.0);

        rSquaredArray.push_back(1 - (ss_res / ss_tot));
    }

    int cFitInd;
    double rSquared;
    for (cFitInd = rSquaredArray.size(); cFitInd >=0; cFitInd-- )
    {
        if (cFitInd == 0)
            return 2;
        if (rSquaredArray[cFitInd] >= cathodicRSquared)
        {
            rSquared = rSquaredArray[cFitInd];
            break;
        }
    }

    std::vector<double> xCathodicLinspace = linspace(xData[0], xData[cFitInd - 1], cFitInd);
    std::vector<double> logXCathodicLinspace = forEachScalar(xCathodicLinspace, 0, [](double a, double b) ->double {return log(a); });
    std::vector<double> yCathodicOpt;
    yCathodicFit.clear();
    cathodicCoeffs.clear();
    std::copy(yData.begin(), yData.begin() + cFitInd, back_inserter(yCathodicOpt));
    polyReg.fitIt(logXCathodicLinspace, yCathodicOpt, 1, cathodicCoeffs);
    yCathodicFit = polynomial(logXCathodicLinspace, cathodicCoeffs);

    double xIntersection = exp((eCorr - cathodicCoeffs[0]) / cathodicCoeffs[1]);

    std::vector<double> anodicCurveCoeffs;
    std::vector<double> xAnodicLinspace = linspace(xAnodic[0], xAnodic[xAnodic.size() - 1], xAnodic.size());
    std::vector<double> logXAnodicLinspace = forEachScalar(xAnodicLinspace, 0, [](double a, double b) ->double {return log(a); });
    polyReg.fitIt(logXAnodicLinspace, yAnodic, 2, anodicCurveCoeffs);
    std::vector<double> yAnodicCurveFit = polynomial(logXAnodicLinspace, anodicCurveCoeffs);

    double step = 0.00001;
    double y2 = eCorr - step;
    std::vector<double> xAnodicFitCoords = {xIntersection, xAnodic[xAnodic.size() - 1]};
    std::vector<double> logXAnodicFitCoords = {log(xAnodicFitCoords[0]), log(xAnodicFitCoords[1])};
    std::vector<double> yAnodicFitCoords = {eCorr, y2};
    double b, a, c;
    int intersects = false;
    std::vector<double> anodicCoeffs;
    std::vector<double> yAnodicFit;
    std::pair<std::complex<double>, std::complex<double>> solutions;
    int realSolutions = false;
    while (!intersects)
    {
        y2 = y2 + step;
        if (y2 >= 0)
            return 3;

        yAnodicFitCoords[1] = y2;
        polyReg.fitIt(logXAnodicFitCoords, yAnodicFitCoords, 1, anodicCoeffs);

        a = anodicCurveCoeffs[2];
        b = anodicCurveCoeffs[1] - anodicCoeffs[1];
        c = anodicCurveCoeffs[0] - anodicCoeffs[0];
        solutions = SolveQuadratic(a, b, c, realSolutions);
        if (realSolutions)
            if (inRange(solutions.first.real(), xAnodicFitCoords[0], xAnodicFitCoords[1]) || inRange(solutions.second.real(), xAnodicFitCoords[0], xAnodicFitCoords[1]))
                intersects = true;
            
    }

    anodicConstant = (yAnodicFitCoords[1] - eCorr) / pow(10, xIntersection - xCathodic[0]);
    cathodicConstant = (eCorr - yCathodicFit[0]) / pow(10, xIntersection - xCathodic[0]);
    iCorr = xAnodicFitCoords[0];
    anodicFitCoords = {iCorr, xAnodicFitCoords[1], eCorr, yAnodicFitCoords[1]};
    cathodicFitCoords = {iCorr, xAnodicFitCoords[1], eCorr, yCathodicFit[0]};

    return 0;
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

    for(int i=0; i < elements - 1; ++i)
    {
      linspaceVector.push_back(startVal + delta * i);
    }
    linspaceVector.push_back(endVal);

    return linspaceVector;
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

