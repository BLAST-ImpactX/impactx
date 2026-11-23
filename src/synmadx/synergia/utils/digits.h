#ifndef DIGITS_H_
#define DIGITS_H_

#include <cmath>

inline long int
digits(long int val)
{
  long int retval = 1;
  long int base = 1;
  for (int i = 0; i < 9; ++i) {
    base *= 10;
    if (val > base) { retval += 1; }
  }
  return retval;
}

inline int
decimal_digits(double val)
{
  int retval = 1;
  double base = 1;
  for (int i = 0; i < 14; ++i) {
    base *= 0.1;
    if (std::fabs(val) < base) { retval += 1; }
  }
  return retval;
}

#endif /* DIGITS_H_ */
