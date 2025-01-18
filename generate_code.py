#!/usr/bin/env python3
import argparse
import math
import os
import shutil
from twiddle_generator import get_nth_root_of_unity_and_psi, twiddle_generator_BR

def is_memory_config_realizable(num_channels, vec_len, coeff_parallel):
    return (num_channels * vec_len) % coeff_parallel == 0

def generate_ini(num_ch, folder):
    if not (1 <= num_ch <= 16):
        raise ValueError("NUM_CH must be between 1 and 16.")

    lines = ["[connectivity]"]

    for i in range(num_ch):
        lines.append(f"sp=ntt_1.x_{i}:HBM[{i}]")

    offset = num_ch

    for i in range(num_ch):
        lines.append(f"sp=ntt_1.y_{i}:HBM[{i + offset}]")

    filename = os.path.join(folder, "link_config.ini")
    with open(filename, "w") as file:
        file.write("\n".join(lines))

    print(f"{filename} has been generated.")

def generate_header(n, mod, bits, B, NUM_CH, folder):
    log2N = int(math.log2(n))
    log2B = int(math.log2(B))

    omega, psi = get_nth_root_of_unity_and_psi(n, mod)

    tw_factors = twiddle_generator_BR(mod, psi, n)

    data_format = None
    if bits == 32:
        data_format = "int"
    else:
        data_format = f"ap_uint<{bits}>"

    # template_file = f"./{folder}/src/ntt.h"
    target_file = os.path.join(folder, "src/ntt.h")

    with open(target_file, "r") as file:
        header_content = file.read()

    header_content = header_content.replace("{DATA_FORMAT}", data_format)
    header_content = header_content.replace("{MOD}", str(mod))
    header_content = header_content.replace("{N}", str(n))
    header_content = header_content.replace("{log2N}", str(log2N))
    header_content = header_content.replace("{B}", str(B))
    header_content = header_content.replace("{log2B}", str(log2B))
    header_content = header_content.replace("{NUM_CH}", str(NUM_CH))
    header_content = header_content.replace("{TW_FACTORS}", ', '.join(map(str, tw_factors)))

    # output_file = os.path.join(folder, "./src/ntt.h")
    with open(target_file, "w") as file:
        file.write(header_content)

    print(f"{target_file} has been generated.")

def generate_makefile(NUM_CH, folder):
    template_file = "./templates/Makefile"
    with open(template_file, "r") as file:
        content = file.read()

    content = content.replace("{NUM_CH}", str(NUM_CH))

    output_file = os.path.join(folder, "Makefile")
    with open(output_file, "w") as file:
        file.write(content)

    print(f"{output_file} has been generated.")

def main():
    parser = argparse.ArgumentParser(description="Calculate the number of NTT cores.")
    parser.add_argument("-N", type=int, default=1024, help="The size of N.")
    parser.add_argument("-q", type=int, default=12289, help="Declare appropriate q")
    parser.add_argument("-bits", type=int, default=32, help="Coefficient bit-width")
    parser.add_argument("-B", type=int, default=8, help="The size of B.")
    parser.add_argument("-NUM_CH", type=int, default=4, help="The number of memory channels (input) ")

    args = parser.parse_args()

    N = args.N
    mod = args.q
    bits = args.bits
    B = args.B
    NUM_CH = args.NUM_CH

    veclen = 512 // bits

    if is_memory_config_realizable(NUM_CH, veclen, 2*B):
        NUM_NTT_CORES = NUM_CH * veclen // (2 * B)

        print(f"Values used -> N: {N}, q: {mod}, bits: {bits}, B: {B}, NUM_CH: {NUM_CH}, veclen: {veclen}")
        print(f"Number of NTT cores: {NUM_NTT_CORES}")

        # Create new folder
        folder_name = f"N_{N}_B_{B}_CH_{NUM_CH}"
        if os.path.exists(folder_name):
            print(f"Folder exists: {folder_name}. No folder created.")
        else:
            print(f"Creating a new folder: {folder_name}")
            os.makedirs(folder_name)

            # Copy all existing files to the new folder
            source_folder = "./templates"
            if os.path.exists(source_folder):
                for item in os.listdir(source_folder):
                    s = os.path.join(source_folder, item)
                    d = os.path.join(folder_name, item)
                    if os.path.isdir(s):
                        shutil.copytree(s, d, dirs_exist_ok=True)
                    else:
                        shutil.copy2(s, d)

            # Generate new files in the folder
            generate_ini(NUM_CH, folder_name)
            generate_header(N, mod, bits, B, NUM_CH, folder_name)
            generate_makefile(NUM_CH, folder_name)

    else:
        print(f"NUM_NTT_CORES is not feasible.")
        return

if __name__ == "__main__":
    main()
