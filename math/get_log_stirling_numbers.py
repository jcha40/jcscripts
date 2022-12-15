import numpy as np


def get_log_stirling_numbers(N):
    """
    Get the log of the unsigned Stirling numbers of the first kind using the recursive formula described here:
        https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind#Recurrence_relation
    Args:
        N: n

    Returns:
        Array of log Stirling numbers for all k given n
    """
    curr = np.zeros(1)
    for n in range(N):
        curr = np.logaddexp(np.log(n) + curr[1:], curr[:-1])  # n * [n k] + [n k-1]
        curr = np.insert(curr, [0, n], [-np.inf, 0])  # put 0 at the beginning and 1 at the end
    return curr
