#!/usr/bin/env python3
import argparse
import math
import os
import shutil
from twiddle_generator import get_nth_root_of_unity_and_psi, twiddle_generator_BR

def check_q_and_data_length(q):
    
    q_list = {12289: 14, 8380417: 23, 3221225473: 32}

    if q not in q_list:
        print(f"Current Q value is not supported")
        return
    
    K = q_list[q]
    
    # Setting up host_data format and bits
    bits = None
    data_format = None

    if K <= 16:
        data_format = "uint16_t"
        bits = 16
    elif K <= 32:
        data_format = "uint32_t"
        bits = 32
    elif K <= 64:
        data_format = "uint64_t"
        bits = 64
    else:
        print(f"Current log q is too large.")
        return
    
    return K, bits, data_format


def round_to_nearest_power_of_2(n):
    if n < 1:
        return 1

    power_low = 2 ** math.floor(math.log2(n))
    power_high = 2 ** math.ceil(math.log2(n))

    return power_low if (n - power_low) < (power_high - n) else power_high


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


def generate_header(n, mod, K, bits, data_format, BU, CH, folder):
    WIDTH = 2*BU
    DEPTH = n / WIDTH
    logDEPTH = int(math.log2(DEPTH))

    log2N = int(math.log2(n))
    log2BU = int(math.log2(BU))

    DATA_BSIZE = int(bits/8)
    DataCHLen = int(64 / DATA_BSIZE)

    DRAM_RATE = 0.5
    EffDataCHLen = round_to_nearest_power_of_2(DataCHLen*DRAM_RATE)

    NUM_CORE = CH * EffDataCHLen / (2*BU)
    GROUP_NUM = int(min(CH, NUM_CORE))
    GROUP_CH_NUM = int(CH / GROUP_NUM)
    MCH = int(GROUP_CH_NUM > 1)

    # it takes too much to calculate psi for large q (3221225473)
    psi = None
    psi_dict = {12289: {64: 140, 128: 8340, 256: 3400, 512: 1987, 1024: 1945}, 
                8380417: {64: 3241972, 128: 1736313, 256: 1921994, 512: 550930, 1024: 1028169},
                3221225473: {64: 1292405718, 128: 1262731197, 256 : 764652596, 512 : 1365964089, 1024: 1168849724} }
    
    if mod in psi_dict:
        if n in psi_dict[mod]:
            psi = psi_dict[mod][n]
        else:
            omega, psi = get_nth_root_of_unity_and_psi(n, mod)

    else:
        print("Current Modulo is not supported: ", mod)
        return

    tw_factors = twiddle_generator_BR(mod, psi, n)

    # template_file = f"./{folder}/src/ntt.h"
    target_file = os.path.join(folder, "src/ntt.h")

    with open(target_file, "r") as file:
        header_content = file.read()
    
    header_content = header_content.replace("{K}", str(K))
    header_content = header_content.replace("{DATA_FORMAT}", data_format)
    header_content = header_content.replace("{MOD}", str(mod))
    header_content = header_content.replace("{N}", str(n))
    header_content = header_content.replace("{log2N}", str(log2N))
    header_content = header_content.replace("{BU}", str(BU))
    header_content = header_content.replace("{log2BU}", str(log2BU))
    header_content = header_content.replace("{logDEPTH}", str(logDEPTH))
    header_content = header_content.replace("{CH}", str(CH))
    header_content = header_content.replace("{MCH}", str(MCH))
    header_content = header_content.replace("{EffDataCHLen}", str(EffDataCHLen))
    header_content = header_content.replace("{DATA_BSIZE}", str(DATA_BSIZE))
    header_content = header_content.replace("{GROUP_NUM}", str(GROUP_NUM))
    header_content = header_content.replace("{GROUP_CH_NUM}", str(GROUP_CH_NUM))
    header_content = header_content.replace("{TW_FACTORS}", ', '.join(map(str, tw_factors)))
    header_content = header_content.replace("{PSI}", str(psi))

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
    # parser.add_argument("-bits", type=int, default=32, help="Coefficient bit-width")
    parser.add_argument("-BU", type=int, default=8, help="The size of BU.")
    parser.add_argument("-CH", type=int, default=4, help="The number of memory channels (input) ")

    args = parser.parse_args()

    N = args.N
    mod = args.q
    # bits = args.bits
    BU = args.BU
    CH = args.CH
    
    K, bits, data_format = check_q_and_data_length(mod)    

    veclen = 512 // bits

    if (CH*veclen) % (2*BU) != 0:
        print(f"{CH} * {veclen} should be divisible by {2 * BU}. Design not feasible.")
        return

    NUM_NTT_CORES = CH * veclen // (2 * BU)

    print(f"Values used -> N: {N}, q: {mod}, HostData: {data_format}, BU: {BU}, CH: {CH}, veclen: {veclen}")
    print(f"Number of NTT cores: {NUM_NTT_CORES}")

    # Create new folder
    folder_name = f"N{N}_BU{BU}_CH{CH}_q{mod}"
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
        generate_header(N, mod, K, bits, data_format, BU, CH, folder_name)
        generate_makefile(CH, folder_name)

 
if __name__ == "__main__":
    main()
