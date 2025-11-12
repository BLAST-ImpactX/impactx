def _kahan_babushka_core(values, return_cumulative=False):
    """Core implementation of the second-order iterative Kahan-Babuska algorithm.

    This is the unified core that implements Klein (2006) algorithm for both
    regular summation and cumulative summation.
    - https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements
    - Klein (2006). "A generalized Kahan–Babuška-Summation-Algorithm". in
      Computing. 76 (3–4). Springer-Verlag: 279–293. doi:10.1007/s00607-005-0139-x

    Args:
        values: Iterable of numeric values to sum
        return_cumulative: If True, returns list of cumulative sums; if False, returns final sum

    Returns:
        float or list: Final sum if return_cumulative=False, list of cumulative sums if True
    """
    sum_val = 0.0
    cs = 0.0  # first-order compensation for lost low-order bits
    ccs = 0.0  # second-order compensation for further lost bits
    c = 0.0  # temporary variable for first-order compensation
    cc = 0.0  # temporary variable for second-order compensation

    if return_cumulative:
        cumulative_sums = [0.0]  # Start with 0.0

    for val in values:
        # First-order Kahan-Babuška step
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

        if return_cumulative:
            # Store the accurate cumulative sum
            cumulative_sums.append(sum_val + cs + ccs)

    if return_cumulative:
        return cumulative_sums
    else:
        return sum_val + cs + ccs


def kahan_babushka_sum(values):
    """Calculate an accurate sum using the second-order iterative Kahan-Babuška algorithm.

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
    return _kahan_babushka_core(values, return_cumulative=False)


def kahan_babushka_cumsum(values):
    """Calculate an accurate cumulative sum using the second-order iterative Kahan-Babuska algorithm.

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
        values: Iterable of numeric values to cumulatively sum

    Returns:
        list: List of cumulative sums with initial 0.0 prepended
    """
    return _kahan_babushka_core(values, return_cumulative=True)
