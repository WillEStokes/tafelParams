#ifndef _SOLVECUBIC_H
#define _SOLVECUBIC_H

#include <utility> // pair
#include <complex> // complex
#include <array>
#include <cmath> // sqrt

int sign(double b);

std::array<std::complex<double>, 3> SolveCubic(double a, double b, double c, double d)
{
  std::array<std::complex<double>, 3> result;

  if (a == 0 || d == 0)
  {
      return result;
  }

  b /= a;
  c /= a;
  d /= a;
  float disc, q, r, dum1, s, t, term1, r13;
  q = (3.0 * c - (b * b)) / 9.0;
  r = -(27.0 * d) + b * (9.0 * c - 2.0 * (b * b));
  r /= 54.0;
  disc = q * q * q + r * r;
  term1 = (b / 3.0);
  if (disc > 0) { // one root real, two are complex
    s = r + sqrt(disc);
    s = ((s < 0) ? -pow(-s, (1.0 / 3.0)) : pow(s, (1.0 / 3.0)));
    t = r - sqrt(disc);
    t = ((t < 0) ? -pow(-t, (1.0 / 3.0)) : pow(t, (1.0 / 3.0)));
    result[0] = -term1 + s + t;
    term1 += (s + t)/2.0;
    result[1] = result[2] = -term1;
    term1 = sqrt(3.0) * (-t + s)/2;
    result[1].imag(term1);
    result[2].imag(-term1);
    return result;
  } 
  
  // The remaining options are all real
  if (disc == 0){ // All roots real, at least two are equal.
    r13 = ((r < 0) ? -pow(-r, (1.0 / 3.0)) : pow(r, (1.0 / 3.0)));
    result[0] = -term1 + 2.0 * r13;
    result[1] = result[2] = -(r13 + term1);
    return result; 
  }
  
  // Only option left is that all roots are real and unequal (to get here, q < 0)
  q = -q;
  dum1 = q * q * q;
  dum1 = acos(r / sqrt(dum1));
  r13 = 2.0 * sqrt(q);

  result[0] = -term1 + r13 * cos(dum1 / 3.0);
  result[1] = -term1 + r13 * cos((dum1 + 2.0 * M_PI) / 3.0);
  result[2] = -term1 + r13 * cos((dum1 + 4.0 * M_PI) / 3.0);

  return result;
}

#endif
