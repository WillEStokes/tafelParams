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
#include "TafelParams.h"

// Function declerations
void readFileWrapper(char* file_name, double data[2817]);
void readFile(std::string file_name, std::vector<double>& data);
std::vector<double> deriv(std::vector<double>& data);
std::vector<double> fastSmooth(std::vector<double>& yData, unsigned int smoothWidth);
std::vector<double> forEachVector(const std::vector<double>& values, const std::vector<double>& modifier, double(*func)(double, double));
std::vector<double> forEachScalar(const std::vector<double>& values, double modifier, double(*func)(double, double));
std::vector<double> polynomial(std::vector<double>& xData, std::vector<double>& coeffs);

int main()
{
    std::vector<double> xDataAll;
    std::vector<double> yDataAll;
    readFile("C:/Users/menwst/Documents/CPP/tafelParams/data/xData", xDataAll);
    readFile("C:/Users/menwst/Documents/CPP/tafelParams/data/yData", yDataAll);

    const double potentialThreshold = 0.01;
    const double cathodicRSquared = 0.9;
    double cathodicConstant;
    double anodicConstant;
    double iCorr;
    double eCorr;
    double cathodicFitCoords[4];
    double anodicFitCoords[4];

    int error = getTafelParams(xDataAll, yDataAll, potentialThreshold, cathodicRSquared, cathodicConstant, anodicConstant, iCorr, eCorr, cathodicFitCoords, anodicFitCoords)
    
}

int getTafelParams(std::vector<double>& xDataAll, 
std::vector<double>& yDataAll, 
const double potentialThreshold, 
const double cathodicRSquared, 
double& cathodicConstant, 
double& anodicConstant, 
double& iCorr, 
double& eCorr, 
double cathodicFitCoords[4], 
double anodicFitCoords[4] )
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
            double eCorr = yData[zeroCrossInd];
            break;
        }
    }
    
    std::vector<double> xCathodic;
    std::vector<double> yCathodic;
    std::vector<double> xAnodic;
    std::vector<double> yAnodic;
    std::copy(xData.begin() + zeroCrossInd + 1, xData.end(), back_inserter(xAnodic));
    std::copy(yData.begin() + zeroCrossInd + 1, yData.end(), back_inserter(yAnodic));

    xCathodic.push_back(xData[0]);
    yCathodic.push_back(yData[0]);
    PolynomialRegression<double> polyReg;
    std::vector<double> logXCathodic;
    std::vector<double> cathodicCoeffs;
    std::vector<double> rSquaredArray;
    std::vector<double> ss_res_array;
    std::vector<double> ss_tot_array;
    std::vector<double> yCathodicFit;
    double yMean, ss_res, ss_tot;
    for (int i = 1; i < zeroCrossInd; i++)
    {
        logXCathodic.clear();
        cathodicCoeffs.clear();
        yCathodicFit.clear();
        
        xCathodic.push_back(xData[i]);
        yCathodic.push_back(yData[i]);

        logXCathodic = forEachScalar(xCathodic, 0, [](double a, double b) ->double {return log(a); });
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

    std::vector<double> xCathodicOpt;
    std::vector<double> yCathodicOpt;
    std::copy(xCathodic.begin(), xCathodic.begin() + cFitInd, back_inserter(xCathodicOpt));
    std::copy(xCathodic.begin(), yCathodic.begin() + cFitInd, back_inserter(yCathodicOpt));
    logXCathodic = forEachScalar(xCathodic, 0, [](double a, double b) ->double {return log(a); });
    polyReg.fitIt(logXCathodic, yCathodic, 1, cathodicCoeffs);
    yCathodicFit = polynomial(logXCathodic, cathodicCoeffs);

    double xIntersection = exp((eCorr - cathodicCoeffs[2]) / cathodicCoeffs[1]);

    std::vector<double> anodicCurveCoeffs;
    std::vector<double> logXAnodic = forEachScalar(xAnodic, 0, [](double a, double b) ->double {return log(a); });
    polyReg.fitIt(logXAnodic, yAnodic, 2, anodicCurveCoeffs);

    double step = 0.001;
    double y2 = eCorr - step;
    std::vector<double> xCoords = {xIntersection, *std::max_element(xAnodic.begin(), xAnodic.end())};
    std::vector<double> yCoords = {eCorr, 0};
    double b, a, c;
    int intersects = false;
    std::vector<double> anodicCoeffs;
    while (!intersects)
    {
        y2 = y2 + step;
        if (y2 >=1)
        {
            return 3;
        }
        yCoords[1] = y2;
        polyReg.fitIt(xCoords, yCoords, 1, anodicCoeffs);
        a = anodicCurveCoeffs[0];
        b = anodicCurveCoeffs[1] + anodicCoeffs[0];
        c = anodicCurveCoeffs[2] + anodicCoeffs[1];
        if (b * b - 4 * a * c >=0)
        {
            intersects = true;
        }
    }

    std::vector<double> yAnodicFit = polynomial(xCoords, anodicCoeffs);

    cathodicConstant = (eCorr - *std::min_element(yCathodicFit.begin(), yCathodicFit.end())) / pow(10, (xIntersection - *std::max_element(xCathodicOpt.begin(), xCathodicOpt.end())));
    anodicConstant = (yAnodicFit[1] - eCorr) / pow(10, (xIntersection - *std::max_element(xCathodicOpt.begin(), xCathodicOpt.end())));
    iCorr = xCoords[0];
    cathodicFitCoords[4] = {iCorr, eCorr, *std::max_element(xCathodicOpt.begin(),  xCathodicOpt.end()), *std::min_element(yCathodicFit.begin(), yCathodicFit.end())};
    anodicFitCoords[4] = {xCoords[0], xCoords[1], yAnodicFit[0], yAnodicFit[1]};

    std::cin.get();
    return 0;
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

