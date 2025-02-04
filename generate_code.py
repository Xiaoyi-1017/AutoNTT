#!/usr/bin/env python3
import argparse
import math
import os
import shutil
from twiddle_generator import get_nth_root_of_unity_and_psi, twiddle_generator_BR

def is_memory_config_realizable(num_channels, vec_len, coeff_parallel):
    return (num_channels * vec_len) % coeff_parallel == 0

def generate_ini(ch, folder):
    if not (1 <= ch <= 16):
        raise ValueError("CH must be between 1 and 16.")

    lines = ["[connectivity]"]
    lines_1 = ["[connectivity]"]

    for i in range(ch):
        lines.append(f"sp=ntt.x_{i}:HBM[{i}]")
        lines_1.append(f"sp=ntt_1.x_{i}:HBM[{i}]")

    offset = ch

    for i in range(ch):
        lines.append(f"sp=ntt.y_{i}:HBM[{i + offset}]")
        lines_1.append(f"sp=ntt_1.y_{i}:HBM[{i + offset}]")

    filename = os.path.join(folder, "link_config.ini")
    with open(filename, "w") as file:
        file.write("\n".join(lines))

    #print(f"{filename} has been generated.")

    filename_1 = os.path.join(folder, "link_config_1.ini")
    with open(filename_1, "w") as file_1:
        file_1.write("\n".join(lines_1))


def generate_header(n, mod, bits, BU, CH, folder):
    WIDTH = 2*BU
    DEPTH = n / WIDTH
    logDEPTH = int(math.log2(DEPTH))

    log2N = int(math.log2(n))
    log2BU = int(math.log2(BU))

    MCH = int((2*BU) > 64 / (bits/8))
    DATA_BSIZE = int(bits/8)
    kVecLen = int(64 / DATA_BSIZE)
    GROUP_NUM = int(min(CH, CH * kVecLen / (2*BU)))
    GROUP_CH_NUM = int(CH / GROUP_NUM)

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
    header_content = header_content.replace("{BU}", str(BU))
    header_content = header_content.replace("{log2BU}", str(log2BU))
    header_content = header_content.replace("{logDEPTH}", str(logDEPTH))
    header_content = header_content.replace("{CH}", str(CH))
    header_content = header_content.replace("{MCH}", str(MCH))
    header_content = header_content.replace("{DATA_BSIZE}", str(DATA_BSIZE))
    header_content = header_content.replace("{GROUP_NUM}", str(GROUP_NUM))
    header_content = header_content.replace("{GROUP_CH_NUM}", str(GROUP_CH_NUM))
    header_content = header_content.replace("{TW_FACTORS}", ', '.join(map(str, tw_factors)))

    # output_file = os.path.join(folder, "./src/ntt.h")
    with open(target_file, "w") as file:
        file.write(header_content)

    #print(f"{target_file} has been generated.")

def generate_makefile(CH, folder):
    template_file = "./templates/Makefile"
    with open(template_file, "r") as file:
        content = file.read()

    content = content.replace("{CH}", str(CH))

    output_file = os.path.join(folder, "Makefile")
    with open(output_file, "w") as file:
        file.write(content)

    #print(f"{output_file} has been generated.")

def main():
    parser = argparse.ArgumentParser(description="Calculate the number of NTT cores.")
    parser.add_argument("-N", type=int, default=1024, help="The size of N.")
    parser.add_argument("-q", type=int, default=12289, help="Declare appropriate q")
    parser.add_argument("-bits", type=int, default=32, help="Coefficient bit-width")
    parser.add_argument("-BU", type=int, default=8, help="The size of BU.")
    parser.add_argument("-CH", type=int, default=4, help="The number of memory channels (input) ")

    args = parser.parse_args()

    N = args.N
    mod = args.q
    bits = args.bits
    BU = args.BU
    CH = args.CH

    veclen = 512 // bits

    if is_memory_config_realizable(CH, veclen, 2*BU):
        NUM_NTT_CORES = CH * veclen // (2 * BU)

        print(f"Values used -> N: {N}, q: {mod}, bits: {bits}, BU: {BU}, CH: {CH}, veclen: {veclen}")
        print(f"Number of NTT cores: {NUM_NTT_CORES}")

        # Create new folder
        folder_name = f"N{N}_BU{BU}_CH{CH}"
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
            generate_ini(CH, folder_name)
            generate_header(N, mod, bits, BU, CH, folder_name)
            generate_makefile(CH, folder_name)

    else:
        print(f"NUM_NTT_CORES is not feasible.")
        return

if __name__ == "__main__":
    main()
