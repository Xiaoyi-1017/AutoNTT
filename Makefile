# Define variables
APP := ntt
HOST := host_ntt

NUM := 100 
TARGET := sw # Default target is 'sw', can be overridden

KERNEL_FILE := ./src/ntt.cpp
HOST_FILE := ./src/host.cpp
CXX_FLAGS := -O2 -fopenmp

PART_NUM := xcu280-fsvh2892-2L-e
CLOCK_PERIOD := 3.33
PLATFORM := xilinx_u280_gen3x16_xdma_1_202211_1

BITSTREAM_XO := $(APP).xo
SIM_DIR1 := ./before_optimized
SIM_DIR2 := ./after_optimized

OPT_WORK_DIR := ./run_optimized

HW_WORK_DIR := ./_x_hw

HW_WORK_DIR_OPT := ./_x_opt

BITSTREAM_HW := $(APP).$(PLATFORM).hw.xclbin

HW_OUTPUT := ./hw_original
HW_OPT_OUTPUT := ./hw_optimized

SAVE_WAVEFORM := 1 # Default to saving waveform

MAKEFLAGS += --no-print-directory

.PHONY: all run clean

# Build target
all:
	@if [ "$(TARGET)" = "sw" ]; then \
		$(MAKE) compile; \
	elif [ "$(TARGET)" = "xo" ]; then \
		$(MAKE) build-xo; \
	elif [ "$(TARGET)" = "opt" ]; then \
		$(MAKE) build-opt; \
	elif [ "$(TARGET)" = "hw" ]; then \
		$(MAKE) build-hw; \
	elif [ "$(TARGET)" = "hw-opt" ]; then \
		$(MAKE) build-hw-opt; \
	else \
		echo "Invalid TARGET: $(TARGET). Use 'sw', 'xo', 'opt', 'hw', or 'hw-opt'."; \
		exit 1; \
	fi

# Run target
run:
	@if [ "$(TARGET)" = "sw" ]; then \
		$(MAKE) run-sw; \
	elif [ "$(TARGET)" = "xo" ]; then \
		$(MAKE) run-xo; \
	elif [ "$(TARGET)" = "opt" ]; then \
		$(MAKE) run-opt; \
	elif [ "$(TARGET)" = "hw" ]; then \
		$(MAKE) run-hw; \
	elif [ "$(TARGET)" = "hw-opt" ]; then \
		$(MAKE) run-hw-opt; \
	else \
		echo "Invalid TARGET: $(TARGET). Use 'sw', 'xo', 'opt', 'hw', or 'hw-opt'."; \
		exit 1; \
	fi

# Compile the application
compile:
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
			-o $(BITSTREAM_XO)
	@echo "XO file generation complete."

# Run XO simulation
run-xo:
	@echo "Running the XO file simulation..."
	@if [ "$(SAVE_WAVEFORM)" -eq 1 ]; then \
		./$(HOST) --bitstream $(BITSTREAM_XO) \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_DIR1) \
			--xosim_save_waveform \
			$(NUM); \
	else \
		./$(HOST) --bitstream $(BITSTREAM_XO) \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_DIR1) \
			$(NUM); \
	fi
	@echo "XO simulation complete."

# Optimize design
build-opt:
	@echo "Optimizing the design..."
	rapidstream gen_config.py
	rapidstream-tapaopt \
		--work-dir $(WORK_DIR) \
		--tapa-xo-path $(BITSTREAM_XO) \
		--device-config ./config_build/device_config.json \
		--floorplan-config ./config_build/floorplan_config.json \
		--connectivity-ini ./link_config.ini
	@echo "Optimization complete."

# Run optimized design
run-opt:
	@echo "Running the optimized design..."
	@if [ "$(SAVE_WAVEFORM)" -eq 1 ]; then \
		./$(HOST) \
			--bitstream "$(WORK_DIR)/dse/solution_0/updated.xo" \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_DIR2) \
			--xosim_save_waveform \
			$(NUM); \
	else \
		./$(HOST) \
			--bitstream "$(WORK_DIR)/dse/solution_0/updated.xo" \
			--xosim-part-num $(PART_NUM) \
			--xosim_work_dir $(SIM_DIR2) \
			$(NUM); \
	fi
	@echo "Optimized design run complete."

# Build hardware configuration
build-hw:
	@echo "Building the hardware configuration..."
	v++ \
		-o $(HW_OUTPUT)/$(BITSTREAM_HW) \
		--temp_dir $(HW_WORK_DIR) \
		--link \
		--target hw \
		--kernel $(APP) \
		--platform $(PLATFORM) \
		--config link_config.ini \
		$(BITSTREAM_XO)
	@echo "Hardware configuration build complete."

build-hw-opt:
	@echo "Building the hardware configuration..."
	v++ \
		-o $(HW_OPT_OUTPUT)/$(BITSTREAM_HW) \
		--temp_dir $(HW_WORK_DIR_OPT) \
		--link \
		--target hw \
		--kernel $(APP) \
		--platform $(PLATFORM) \
		--config link_config.ini \
		$(WORK_DIR)/dse/solution_0/updated.xo
	@echo "Hardware configuration build complete."

# Run on hardware
run-hw:
	@echo "Running on hardware..."
	./$(HOST) --bitstream $(HW_OUTPUT)/$(BITSTREAM_HW) $(NUM)
	@echo "Hardware run complete."

run-hw-opt:
	@echo "Running on hardware..."
	./$(HOST) --bitstream $(HW_OPT_OUTPUT)/$(BITSTREAM_HW) $(NUM)
	@echo "Hardware run complete."

# Clean generated files
clean:
	@echo "Cleaning up..."
	@rm -rf $(HOST) $(BITSTREAM_XO)
	@rm -rf work.out $(OPT_WORK_DIR) $(HW_WORK_DIR) $(HW_OUTPUT) $(HW_WORK_DIR_OPT) $(HW_OPT_OUTPUT)
	@rm -rf $(SIM_DIR1) $(SIM_DIR2)