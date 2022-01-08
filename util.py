from random import randrange

def miller_rabin(n, a = None):
    """
    Perform the Miller-Rabin primality testing algorithm.

    Returns False if n is determined to be composite.
    Returns True if n is a probable prime.

    Taken and ever so slightly modified from 
    https://github.com/jaanos/kirv/blob/master/algorithms/primality.py
    """
    assert n > 1
    m = n-1
    k = 0
    while m % 2 == 0:
        m //= 2
        k += 1

    if a is None:
        a = randrange(2, n-1)

    b = pow(a, m, n)
    if b % n == 1:
        return True

    for _ in range(k):
        if (b+1) % n == 0:
            return True
        b = b*b % n

    return False # composite

def inverse(a, n, strict = True):
    """
    Compute the inverse of a modulo n using the Extended Euclidean algorithm.
    If strict is True (default), raises an error if a is not invertible.
    Otherwise, a number b is output such that a*b % n == gcd(a, b).
    """
    b, x, y = n, 0, 1
    while a != 0:
        k = int(b // a)
        b, a, x, y = a, b - k*a, y, x - k*y
    if strict and b not in [-1, 1]:
        raise ValueError("input not invertible")
    return (b*x) % n

def bit_length(x):
    """
    Number of bits in the binary representation of x.
    """
    return len(bin(x)) - 2