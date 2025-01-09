#!/usr/bin/env python3
import argparse
import math
import os
from twiddle_generator import get_nth_root_of_unity_and_psi, twiddle_generator_BR

def is_memory_config_realizable(num_channels, vec_len, coeff_parallel):
    
    return (num_channels * vec_len) % coeff_parallel == 0

def generate_ini(num_ch, filename="link_config.ini"):

    if not (1 <= num_ch <= 16):
        raise ValueError("NUM_CH must be between 1 and 16.")

    lines = ["[connectivity]"]

    for i in range(num_ch):
        lines.append(f"sp=ntt_1.x_{i}:HBM[{i}]")

    offset = num_ch

    for i in range(num_ch):
        lines.append(f"sp=ntt_1.y_{i}:HBM[{i + offset}]")

    with open(filename, "w") as file:
        file.write("\n".join(lines))

    print(f"{filename} has been generated.")



def generate_header(n, mod, bits, B, NUM_CH, output_file="./src/ntt.h"):
    # Calculating necessary numbers
    log2N = int(math.log2(n))
    log2B = int(math.log2(B))

    omega, psi = get_nth_root_of_unity_and_psi(n, mod)

    tw_factors = twiddle_generator_BR(mod, psi, n)

    data_format = None
    if bits == 32:
        data_format = "int"
    else:
        data_format = f"ap_uint<{bits}>"

    # Read the template file
    template_file = "./templates/ntt.h"
    with open(template_file, "r") as file:
        header_content = file.read()

    # Replace placeholders in the template with actual values
    header_content = header_content.replace("{DATA_FORMAT}", data_format)
    header_content = header_content.replace("{MOD}", str(mod))
    header_content = header_content.replace("{N}", str(n))
    header_content = header_content.replace("{log2N}", str(log2N))
    header_content = header_content.replace("{B}", str(B))
    header_content = header_content.replace("{log2B}", str(log2B))
    header_content = header_content.replace("{NUM_CH}", str(NUM_CH))
    header_content = header_content.replace("{TW_FACTORS}", ', '.join(map(str, tw_factors)))


    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)


    # Write the modified header to the output file
    with open(output_file, "w") as file:
        file.write(header_content)

    print(f"{output_file} has been generated.")


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Calculate the number of NTT cores.")
    parser.add_argument("-N", type=int, default=1024, help="The size of N.")
    parser.add_argument("-q", type=int, default=12289, help="Declare appropriate q")
    parser.add_argument("-bits", type=int, default=32, help="Coefficient bit-width")
    parser.add_argument("-B", type=int, default=8, help="The size of B.")
    parser.add_argument("-NUM_CH", type=int, default=4, help="The number of memory channels (input) ")

    args = parser.parse_args()

    # Assign arguments to variables
    N = args.N
    mod = args.q
    bits = args.bits
    B = args.B
    NUM_CH = args.NUM_CH
    
    # Calculate vector length
    veclen = 512 // bits


    # Check if NUM_NTT_CORES is an integer
    if is_memory_config_realizable(NUM_CH, veclen, 2*B):

        NUM_NTT_CORES = NUM_CH * veclen // (2 * B)

        print(f"Values used -> N: {N}, B: {B}, NUM_CH: {NUM_CH}, veclen: {veclen}")
        print(f"Number of NTT cores: {NUM_NTT_CORES} ")
        print(f"number of channels per ntt_core: {(2*B) // veclen}")
        print(f"number of ntt_cores per channel: {veclen // (2*B)}")

        # Generate link_config.ini
        generate_ini(NUM_CH)
        
        # Generate ntt.h
        generate_header(N, mod, bits, B, NUM_CH)

    
    else:
        print(f"NUM_NTT_CORES is not feasible: {NUM_NTT_CORES}")
        return

    

if __name__ == "__main__":
    main()
