#ifndef _QUADRATICSOLUTIONS_H
#define _QUADRATICSOLUTIONS_H

#include <utility> // pair
#include <complex> // complex
#include <cmath> // sqrt

int sign(double b);

std::pair<std::complex<float>, std::complex<float>> SolveQuadratic(double a, double b, double c)
{
  std::pair<std::complex<float>, std::complex<float>> result;

  if(a < 0.000001) // ==0
  {
    result.first = result.second = NAN;
    return result;
  }

  double delta = b * b - 4 * a * c;
  if(delta >= 0)
  {
    result.first = (-0.5 * (b + sign(b) * sqrt(delta))) / a;
    result.second = c / (-0.5 * (b + sign(b) * sqrt(delta)));
  }
  else
  {
    result.first = result.second = -b / 2 / a;
    result.first.imag(sqrt(-delta) / 2 / a);
    result.second.imag(-sqrt(-delta) / 2 / a);
  }

  return result;
}

int sign(double b)
{
  if (b >= 0)
    return 1;
  else
    if (b < 0)
    return -1;
  return 0;
}

#endif
