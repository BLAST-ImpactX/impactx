#ifndef FAST_INT_FLOOR_H_
#define FAST_INT_FLOOR_H_

inline int
fast_int_floor(const double x)
{
  int ix = static_cast<int>(x);
  return x > 0.0 ? ix : ((x - ix == 0) ? ix : ix - 1);
}

#endif /* FAST_INT_FLOOR_H_ */
