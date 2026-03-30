#ifndef PERIOD_H_
#define PERIOD_H_

template <class T>
class bunch_t;

using Bunch = bunch_t<double>;

void apply_longitudinal_bucket_barrier(Bunch& bunch, double length);
void apply_longitudinal_periodicity(Bunch& bunch, double length);
void apply_zcut(Bunch& bunch, double length);

#endif /* PERIOD_H_ */
