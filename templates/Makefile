# Define variables
APP := ntt
HOST := host_ntt

NUM := 100 

# Default TARGET is 'sw', and it can be overridden by specifying a different value
TARGET := sw

KERNEL_FILE := ./src/ntt.cpp
HOST_FILE := ./src/host.cpp
CXX_FLAGS := -O2 -fopenmp

PART_NUM := xcu280-fsvh2892-2L-e
CLOCK_PERIOD := 3.33
PLATFORM := xilinx_u280_gen3x16_xdma_1_202211_1

SIM_XO_DIR := ./sim_xo
SIM_XO_OPT_DIR := ./sim_xo_opt
OPT_DIR := ./opt
CONFIG_BUILD_DIR := $(OPT_DIR)/config_build
OPT_WORK_DIR := ./x_opt
HW_WORK_DIR := ./x_hw

XO_FILE := $(APP).xo
BITSTREAM_HW := $(APP).$(PLATFORM).hw.xclbin
BITSTREAM_HW_OPT := $(APP)_$(PLATFORM).xclbin

HW_OUTPUT := ./hw
HW_OPT_OUTPUT := ./hw_opt

# Default is to save waveform
SAVE_WAVEFORM := 1 

MAKEFLAGS += --no-print-directory

NUM_CH := {NUM_CH} 

.PHONY: all compile run clean help

VALID_TARGETS := sw xo xo-opt hw hw-opt

# Build target
all:
	@if echo "$(VALID_TARGETS)" | grep -wq "$(TARGET)"; then \
		$(MAKE) build-$(TARGET); \
	else \
		echo "Invalid TARGET: $(TARGET). Use one of: $(VALID_TARGETS)."; \
		exit 1; \
	fi

# Run target
run:
	@if echo "$(VALID_TARGETS)" | grep -wq "$(TARGET)"; then \
		$(MAKE) run-$(TARGET); \
	else \
		echo "Invalid TARGET: $(TARGET). Use one of: $(VALID_TARGETS)."; \
		exit 1; \
	fi

# Compile the application
build-sw compile:
	@echo "Compiling the application..."
	tapa g++ $(KERNEL_FILE) $(HOST_FILE) -o $(HOST) $(CXX_FLAGS)
	@if [ $$? -eq 0 ]; then \
		echo "Compile successful!"; \
	else \
		echo "Compile failed!" >&2; \
		exit 1; \
	fi

# Run the compiled application @echo "Running the compiled application..."
run-sw:
	@if [ ! -f $(HOST) ]; then \
		echo "Error: Executable $(HOST) not found. Please build it first!" >&2; \
		exit 1; \
	fi
	./$(HOST) $(NUM)

# Build XO file
build-xo:
	@echo "Building the XO file..."
		tapa compile \
			--top $(APP) \
			--part-num $(PART_NUM) \
			--clock-period $(CLOCK_PERIOD) \
			-f $(KERNEL_FILE) \
			-o $(XO_FILE)
	@echo "XO file generation complete."

# Run XO simulation
run-xo:
	@echo "Running the XO file simulation..."
	@if [ "$(SAVE_WAVEFORM)" -eq 1 ]; then \
		./$(HOST) --bitstream $(XO_FILE) \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_XO_DIR) \
			--xosim_save_waveform \
			$(NUM); \
	else \
		./$(HOST) --bitstream $(XO_FILE) \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_XO_DIR) \
			$(NUM); \
	fi
	@echo "XO simulation complete."

# Optimize design and build xo file
build-xo-opt:
	@echo "Making configuration files for rapidstream..."
	mkdir -p $(OPT_DIR)
	(cd $(OPT_DIR) && rapidstream ../gen_config.py --num_ch $(NUM_CH))
	@echo "Optimizing the design..."
	rapidstream-tapaopt \
		--work-dir $(OPT_DIR) \
		--tapa-xo-path $(XO_FILE) \
		--device-config $(CONFIG_BUILD_DIR)/device_config.json \
		--floorplan-config $(CONFIG_BUILD_DIR)/floorplan_config.json \
		--connectivity-ini ./link_config_1.ini
	@echo "Optimization complete."

# Optimize design and build FPGA bitstream
build-hw-opt:
	@echo "Making configuration files for rapidstream..."
	mkdir -p $(OPT_DIR)
	(cd $(OPT_DIR) && rapidstream ../gen_config.py --num_ch $(NUM_CH))
	@echo "Optimizing the design..."
	rapidstream-tapaopt \
		--work-dir $(OPT_DIR) \
		--tapa-xo-path $(XO_FILE) \
		--device-config $(CONFIG_BUILD_DIR)/device_config.json \
		--floorplan-config $(CONFIG_BUILD_DIR)/floorplan_config.json \
		--implementation-config ./impl_config.json \
		--connectivity-ini ./link_config.ini \
		--run-impl
	@echo "Optimization complete."

# Run optimized design
run-xo-opt:
	@echo "Running the optimized design..."
	@if [ "$(SAVE_WAVEFORM)" -eq 1 ]; then \
		./$(HOST) \
			--bitstream "$(OPT_DIR)/dse/solution_0/updated.xo" \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_XO_OPT_DIR) \
			--xosim_save_waveform \
			$(NUM); \
	else \
		./$(HOST) \
			--bitstream "$(OPT_DIR)/dse/solution_0/updated.xo" \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_XO_OPT_DIR) \
			$(NUM); \
	fi
	@echo "Optimized design run complete."

# Build the FPGA bitstream
build-hw:
	@echo "Building the hardware configuration..."
	v++ \
		-o $(HW_OUTPUT)/$(BITSTREAM_HW) \
		--temp_dir $(HW_WORK_DIR) \
		--link \
		--target hw \
		--kernel $(APP) \
		--platform $(PLATFORM) \
		--config link_config_1.ini \
		$(XO_FILE)
	@echo "Hardware configuration build complete."

build-hw-opt-old:
	@echo "Building the hardware configuration..."
	v++ \
		-o $(HW_OPT_OUTPUT)/$(BITSTREAM_HW) \
		--temp_dir $(OPT_WORK_DIR) \
		--link \
		--target hw \
		--kernel $(APP) \
		--platform $(PLATFORM) \
		--config link_config.ini \
		$(OPT_DIR)/dse/solution_0/updated.xo
	@echo "Hardware configuration build complete."

# Run on FPGA
run-hw:
	@echo "Running on hardware..."
	./$(HOST) --bitstream $(HW_OUTPUT)/$(BITSTREAM_HW) $(NUM)
	@echo "Hardware run complete."

run-hw-opt:
	@echo "Running on hardware..."
	./$(HOST) --bitstream $(OPT_DIR)/dse/solution_0/vitis_run_hw/$(BITSTREAM_HW_OPT) $(NUM)
	@echo "Hardware run complete."

run-hw-opt-old:
	@echo "Running on hardware..."
	./$(HOST) --bitstream $(HW_OPT_OUTPUT)/$(BITSTREAM_HW) $(NUM)
	@echo "Hardware run complete."

#Help
help:
	@echo "Usage:"
	@echo "  Build Commands:"
	@echo "    make all TARGET=<type>  # Build the specified TARGET"
	@echo "      Valid TARGET values:"
	@echo "        sw      - Build the host executable"
	@echo "        xo      - Build the kernel XO file"
	@echo "        hw      - Build the FPGA bitstream"
	@echo "        hw-opt  - Optimize with Rapidstream and build FPGA bitstream"
	@echo ""
	@echo "  Run Commands:"
	@echo "    make run TARGET=<type> <sample_num>  # Run with the TARGET for specified number of samples"
	@echo "      Valid TARGET values:"
	@echo "        sw      - Run software simulation"
	@echo "        xo      - Simulate with the XO file"
	@echo "        hw      - Run the FPGA bitstream"
	@echo "        hw-opt  - Run the RS-optimized FPGA bitstream"
	@echo ""

# Clean generated files
clean:
	@echo "Cleaning up..."
	@rm -f $(HOST) $(XO_FILE)
	@rm -rf work.out .ipcache .Xil $(OPT_DIR) $(HW_WORK_DIR) $(HW_OUTPUT) $(OPT_WORK_DIR) $(HW_OPT_OUTPUT)
	@rm -rf $(SIM_XO_DIR) $(SIM_XO_OPT_DIR)
	@rm -f v++*.log xcd.log gurobi.log
	@echo "Clean up complete."
