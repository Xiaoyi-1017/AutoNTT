# AutoNTT
Code generator for FPGA NTT accelerator

### Prerequisites

- AMD Vitis (2023.2) - https://www.amd.com/en/products/software/adaptive-socs-and-fpgas/vitis.html

- AMD Alveo U280 Platform - https://www.xilinx.com/publications/product-briefs/alveo-u280-product-brief.pdf

- TAPA (0.1.20250329) - https://tapa.readthedocs.io/en/main/

- Rapidstream (2025.1.0212, recommended) - https://docs.rapidstream-da.com/

### Install Python Requirements
```bash
pip install -r requirements.txt 
```

## Generating an NTT project 

Please use generate_code.py to automatically generate an NTT project.

| Argument | Description | Supported Values | Default |
|----------|-------------|-----------------|---------|
| **N** | Transform size (polynomial degree) | 256, 512, 1024 | 1024 |
| **q** | Prime modulus | 12289, 8380417, 3221225473 | 12289 |
| **BU** | Number of butterfly units per stage | 1, 2, 4, 8, 16, 32 | 16 |
| **CH** | Number of input HBM channels <br> (Total of 2×CH channels are used for input & output) | 1, 2, 4, 8 | 8 |

Example command (with default values):
```bash
./generate_code.py -N 1024 -q 12289 -BU 16 -CH 8 
```

Example console message:
```
Values used -> N: 1024, q: 12289, HostData: uint16_t, BU: 16, CH: 8, veclen: 32
Number of NTT cores: 8
Creating a new folder: N1024_BU16_CH8_q12289
```

## Compilation & Execution 

You can use the Makefile for compilation & execution.
We recommend that you use Rapidstream to achieve higher clock frequency, but you could only use the traditional Vitis flow.

HOST Compilation
```bash
make all TARGET=sw
```
Synthesize and generate xo file
```bash
make all TARGET=xo
```
Generate bitstream from xo file (Vitis only, no Rapidstream)
```bash
make all TARGET=hw
```
Optimize xo file with Rapidstream and generate bitstream (recommended)
```bash
make all TARGET=hw-opt
```
Simulation / Onboard verification (NUM = number of polynomials)
```bash
make run TARGET={sw | xo | hw | hw-opt} NUM=10000
```

You should see **PASSED!** message after running on-board verification.

