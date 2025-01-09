# NTT_accel_generator
Code generator for NTT accelerator

### Prerequisites

AMD Vitis (2022.2) - https://www.amd.com/en/products/software/adaptive-socs-and-fpgas/vitis.html

TAPA - https://tapa.readthedocs.io/en/main/

Rapidstream (Not used) - https://docs.rapidstream-da.com/

### Install Python Requirements
```bash
pip install -r requirements.txt 
```

## How to use MAKEFILE
Actions (TARGET = sw, xo, hw)
1) all: Host compilation [tapa g++] / Synthesis into RTL [tapa compile] / Bitstream generation [v++]
2) run: sw simulation / RTL simulation / onboard


HOST Compilation
```bash
make all TARGET=sw
```
Synthesis into RTL
```bash
make all TARGET=xo
```

Bitstream generation
```bash
make all TARGET=hw
```

Simulation / Onboard verification (NUM = number of samples)
```bash
make run TARGET=(sw / xo / hw)  NUM=10000
```

## How to use generate_code.py

Checking if B (number of butterfly units per NTT CORE) and NUM_CH (Number of memory channels) are compatible

Condition: NUM_CH * VECLEN / number of coefficients in parallel -> Integer

Arguments
- **N**: transform size, polynomial degree
- **q**: prime modulus (i.e. 12289)
- **bits**: bit length of coefficients
- **B**: number of butterfly units in parallel (per stage)
- **NUM_CH**: number of memory banks used for input



FULL CODE (DEFAULT)
```bash
./generate_code.py -N 1024 -q 12289 -bits 32 -B 8 -NUM_CH 4 
```
Example 1: B = 16, NUM_CH = 2
```bash
./generate_code.py -B 16 -NUM_CH 2 
```
Terminal output
```bash
Values used -> N: 1024, B: 8, NUM_CH: 4, veclen: 16
Number of NTT cores: 4 
number of channels per ntt_core: 1
link_config.ini has been generated.
./src/ntt.h has been generated.
./src/ntt.cpp has been generated. 
```


## Updates - 2025/01/10
### Modified
- **ntt.cpp**: implemented input_stage and temporal_stage differently so that it does not require further modification
- **host.cpp**: implemented input & output processing that supports B less than 8 (multiple samples in the same channel)
- **generate_code.py**: creates folder containing necessary files (src, Makefile, link_config.ini) for different configurations 

### Added
- **gen_config.py**: generates floorplan & device configs for Rapidstream-tapaopt (stored under "all_files" folder)
