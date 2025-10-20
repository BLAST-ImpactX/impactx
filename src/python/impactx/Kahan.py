def kahan_babushka_sum(values):
    """Calculate an accurate sum using the second-order iterative Kahan-Babuska algorithm.

    This implementation follows Klein (2006) to provide high-precision summation
    that avoids floating-point precision errors when summing many small values.
    - https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements
    - Klein (2006). "A generalized Kahan–Babuška-Summation-Algorithm". in
      Computing. 76 (3–4). Springer-Verlag: 279–293. doi:10.1007/s00607-005-0139-x

    The algorithm uses second-order compensation for lost low-order bits during
    floating-point addition, providing significantly better accuracy than naive
    summation when dealing with large numbers of small values (e.g., many ds
    values in a long lattice).

    Args:
        values: Iterable of numeric values to sum

    Returns:
        float: Accurate sum of all values
    """
    sum_val = 0.0
    cs = 0.0  # first-order compensation for lost low-order bits
    ccs = 0.0  # second-order compensation for further lost bits

    for val in values:
        # First-order Kahan-Babuska step
        t = sum_val + val
        if abs(sum_val) >= abs(val):
            c = (sum_val - t) + val
        else:
            c = (val - t) + sum_val
        sum_val = t

        # Second-order compensation step
        t = cs + c
        if abs(cs) >= abs(c):
            cc = (cs - t) + c
        else:
            cc = (c - t) + cs
        cs = t
        ccs += cc

    return sum_val + cs + ccs
