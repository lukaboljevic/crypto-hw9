
from math import ceil
from random import randrange
from hashlib import sha256
from util import miller_rabin, inverse, bit_length

acceptable_pairs = [(1024, 160), (2048, 224), (2048, 256), (3072, 256)]
pair_to_iteration = {
    (1024, 160): 40,
    (2048, 224): 56,
    (2048, 256): 56,
    (3072, 256): 64
}
# Taken from https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf
#     Part C.3, Table C.1


def hash_and_convert(input: str):
    """
    Perform SHA256 on the input and return the decimal form of the result
    """
    digest = sha256(input.encode()).hexdigest()
    return int(digest, 16)

def prime(number, L, N):
    """
    Perform the Miller Rabin test enough times for the given pair L, N. 

    Return False if Miller Rabin fails at least once for the given number
    (i.e. says that the number is composite)
    Return True if for all of the iterations, Miller Rabin determined
    this number is prime
    """
    num_iterations = pair_to_iteration[(L, N)]
    for _ in range(num_iterations):
        if not miller_rabin(number):
            return False
    return True

def generate_probable_primes(L, N, seedlen):
    """
    Generate two (very highly) probable primes p and q, of bit length L and N respectively.

    Taken from https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf
        Part A.1.1.2

    Parameters
    ----------
    L: int
        Desired length of prime p (in bits).
    N: int
        Desired length of prime q (in bits).
    seedlen: int
        Desired length of the domain_parameter_seed.
        Must be greater than or equal to N.

    Returns
    -------
    p: int
        -1 if generation failed.
    q: int
        -1 if generation failed.
    domain_parameter_seed: int
    """
    invalid = -1, -1, -1

    # Check simple conditions, steps 1 and 2
    if (L, N) not in acceptable_pairs: 
        return invalid
    if seedlen < N: 
        return invalid
    
    # Generate once all the necessary stuff
    # Steps 3, 4
    outlen = seedlen*2 # for example
    n = ceil(L / outlen) - 1
    b = L - 1 - (n * outlen)
    seed_lower_bound = pow(2, seedlen - 1)
    seed_upper_bound = pow(2, seedlen)
    two_to_N_minus_1 = pow(2, N-1)
    two_to_outlen = pow(2, outlen)
    two_to_seedlen = pow(2, seedlen)
    two_to_b = pow(2, b)
    two_to_L_minus_1 = pow(2, L-1)
    visited = set()

    # Try to generate p and q
    while True:
        # Step 5 technically
        print("Attempting to generate p and q...")
        print("Generating q...")
        while True:
            domain_parameter_seed = randrange(seed_lower_bound, seed_upper_bound)
            if domain_parameter_seed in visited:
                # don't try the same ones
                continue
            visited.add(domain_parameter_seed)

            # Step 6
            decimal_hash = hash_and_convert(str(domain_parameter_seed))
            U = decimal_hash % two_to_N_minus_1

            # Step 7
            q = two_to_N_minus_1 + U + 1 - (U % 2)
            
            # Step 8, 9
            if prime(q, L, N):
                break

        print("q successfully generated!")
        print("Generating p...")
        visited.clear() # clear for later, if we fail to find p

        # Step 10
        offset = 1

        # Step 11
        for counter in range(4*L): # [0, 4L-1]

            # Step 11.1 and 11.2 at the same time 
            W = 0
            powers_two = two_to_outlen # 2^outlen
            for j in range(n+1): # [0, n]
                temp = domain_parameter_seed + offset + j
                temp %= two_to_seedlen
                v_j = hash_and_convert(str(temp))
                if j == n:
                    v_j %= two_to_b

                if j >= 2:
                    powers_two *= 2 # (2^outlen)^2, (2^outlen)^3, (2^outlen)^4, ...

                if j >= 1:
                    v_j *= powers_two # v_j *= (2^outlen)^j

                W += v_j
            
            # Step 11.3
            X = W + two_to_L_minus_1

            # Step 11.4
            c = X % (2*q)

            # Step 11.5
            p = X - (c - 1)

            # Step 11.6, 11.7, 11.8
            if p >= two_to_L_minus_1:
                if prime(p, L, N):
                    print("Yay!! Successfully generated p and q")
                    print()
                    return p, q, domain_parameter_seed

            # Step 11.9
            offset = offset + n + 1
        
        # Step 12 more or less
        print("Failed to generate p =(")
        print("Retrying...")
        print("================================================")

def generate_generator(p, q, domain_parameter_seed, index):
    """
    Generate a generator g of order q, in the multiplicative group Z_p*.

    Taken from https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf
        Part A.2.3

    Parameters:
    -----------
    p: int
        Prime p.
    q: int
        Prime q.
    domain_parameter_seed: int
        The domain parameter seed used during the generation of p and q.
    index: int
        The index to be used for generating g. This is a bit string of
        length 8 that represents an unsigned integer. (So, integer from [0, 255])

    Returns
    -------
    g: int
        The generated generator.
    """
    # count is an unsigned 16 bit integer i.e. integer in range [0, 65,535]

    # Step 1
    if index < 0 or index > 255:
        return -1

    # Step 3
    e = (p - 1) // q # q is a prime divisor of p-1

    # Convert all the stuff needed for step 7 to binary
    ggen_dec = int("0x6767656E", 16)
    ggen_bin = bin(ggen_dec)[2:] # eliminate 0b
    domain_parameter_seed_bin = bin(domain_parameter_seed)[2:]
    index_bin = bin(index)[2:]

    # Step 4
    count = 0

    # print("Attempting to generate g...")
    # Step 5
    while True:

        # Step 6
        if count == 65535: # (2^16 - 1)
            # because increasing by 1 will (would) wrap it to 0
            print("Failed to generate g =(")
            print("Try again")
            print("================================================")
            return -1
        count += 1
        count_bin = bin(count)[2:]

        # Step 7
        U = domain_parameter_seed_bin + ggen_bin + index_bin + count_bin

        # Step 8
        W = hash_and_convert(U)

        # Step 9
        g = pow(W, e, p)

        # Steps 10 & 11
        if g >= 2:
            # print("Yay!! Successfully generated g")
            # print()
            return g

def generate_k(p, q, g):
    """
    Generate the secret number k (and its inverse) unique to each message. 
    This k and its inverse are then used for signing a message.

    Taken and modified from
    https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf
        Part B.2.2

    Parameters:
    -----------
    p: int
    q: int
    g: int
        The two primes and generator element.

    Returns:
    --------
    k: int
        Secret number k used for signing
    k_inv: int
        The inverse of k
    """
    invalid = (-1, -1)
    L = bit_length(p)
    N = bit_length(q)
    print(L, N)
    if (L, N) not in acceptable_pairs:
        return invalid

    # https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-57pt1r5.pdf
    # page 67

    lower_bound = pow(2, N//2 - 1)
    upper_bound = (q - 2) + 1
    visited = set()

    while True:    
        c = randrange(lower_bound, upper_bound) # has at least half as many bits as N, not bigger than q-2 (that's how I understood it)
        if c in visited:
            continue
        visited.add(c)
        k = c + 1 # not bigger than q-1
        lhs = pow(g, k, q)
        rhs = pow(g, k, p) % q
        if lhs == rhs:
            print("WHOOPS")
            continue
        k_inv = inverse(k, q)
        return k, k_inv


L = 1024
N = 160
seedlen = 170
p, q, dps = generate_probable_primes(L, N, seedlen)
assert p != -1 and q != -1 and dps != -1

g = generate_generator(p, q, dps, 1)
assert g != -1

k, k_inv = generate_k(p, q, g)
assert k != -1 and k_inv != -1

print(f"p: {p}\n")
print(f"q: {q}\n")
print(f"g: {g}\n")
print(f"k: {k}\n")
print(f"k_inv: {k_inv}")