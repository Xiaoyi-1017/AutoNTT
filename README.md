# NTT_accel_generator
Code generator for NTT accelerator

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