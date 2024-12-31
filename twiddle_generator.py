import numpy as np
from sympy import primitive_root, mod_inverse


def reverse_bits(num, bits):
    """
    Reverse the bits of an integer `num` within the range defined by `bits`.
    """
    reversed_num = 0
    for i in range(bits):
        if num & (1 << i):
            reversed_num |= 1 << ((bits - 1) - i)
    return reversed_num

def bit_reverse_array(arr):
    """
    Perform bit-reversal ordering on the input array `arr`.
    """
    n = len(arr)
    bits = int(np.log2(n))
    out = np.zeros_like(arr)
    for i in range(n):
        out[reverse_bits(i, bits)] = arr[i]
    return out


def get_nth_root_of_unity(n, mod):
    """
    Calculate the nth root of unity (omega) modulo `mod` using sympy.
    """
    g = primitive_root(mod)  # Find a primitive root of the modulus
    omega = pow(g, (mod - 1) // n, mod)  # Calculate the nth root of unity
    return omega


def get_psi(n, q, omega):
    """
    Find a suitable psi value that satisfies the conditions.
    """
    for psi in range(q):
        mul = pow(psi, 2 * n, q)
        if mul == 1 and pow(psi, 2, q) == omega % q:
            return psi
    return None


def get_nth_root_of_unity_and_psi(n, mod):
    """
    Calculate the nth root of unity (omega) modulo `mod` using sympy.
    """
    
    g = primitive_root(mod)  # Find a primitive root of the modulus
    omega = pow(g, (mod - 1) // n, mod)  # Calculate the nth root of unity


    psi = get_psi(n, mod, omega)
    
    return omega, psi



def twiddle_generator_BR(mod, psi, n):
    """
    Generate bit-reversed twiddle factors.
    """
    arr = [1]*n

    for i in range(n-1):
        arr[i+1] = (arr[i]*psi) % mod

    return bit_reverse_array(arr)



def main():
    n = 2048
    mod = 12289
    omega = get_nth_root_of_unity(n, mod)
    psi = get_psi(n, mod, omega)    

    
    print("Omega:", omega)
    print("Psi:", psi)
    
    # Generate bit-reversed twiddle factors
    tw_factors_BR = twiddle_generator_BR(mod, psi, n)

    print("twiddle factors: \n", tw_factors_BR)


if __name__ == "__main__":
    main()
