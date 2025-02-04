"""Generate rapidstream configuration files."""
__copyright__ = """
Copyright (c) 2024 RapidStream Design Automation, Inc. and contributors.
All rights reserved. The contributor(s) of this file has/have agreed to the
RapidStream Contributor License Agreement.
"""
import os
from pathlib import Path
import argparse

#from rapidstream.assets.device_library.u280.u280 import get_u280_precollected_device
from rapidstream import get_u280_vitis_device_factory, RapidStreamTAPA
from rapidstream.assets.floorplan.floorplan_config import FloorplanConfig
from rapidstream import PipelineConfig

CONFIG_DIR = Path("config_build")
FLOOR_PLAN_CONFIG = CONFIG_DIR / "floorplan_config.json"
DEVICE_CONFIG = CONFIG_DIR / "device_config.json"
PIPELINE_CONFIG = CONFIG_DIR / "pipeline_config.json"
VITIS_PLATFORM = "xilinx_u280_gen3x16_xdma_1_202211_1"

def gen_config(ch: int) -> None:
    """Generate configuration files."""
    if not os.path.exists(CONFIG_DIR):
        os.makedirs(CONFIG_DIR)

    port_pre_assignments = {
        "ap_clk": "SLOT_X0Y0:SLOT_X0Y0",
        "ap_rst_n": "SLOT_X0Y0:SLOT_X0Y0",
        "interrupt": "SLOT_X0Y0:SLOT_X0Y0",
        "s_axi_control_.*": "SLOT_X0Y0:SLOT_X0Y0",
    }
    for i in range(ch):
        if i < 8:
            port_pre_assignments[f".*m_axi_x_{i}_.*"] = "SLOT_X0Y0:SLOT_X0Y0"
            port_pre_assignments[f".*m_axi_y_{i}_.*"] = "SLOT_X0Y0:SLOT_X0Y0"
        else:
            port_pre_assignments[f".*m_axi_x_{i}_.*"] = "SLOT_X1Y0:SLOT_X1Y0"
            port_pre_assignments[f".*m_axi_y_{i}_.*"] = "SLOT_X1Y0:SLOT_X1Y0"


    floorplan_config = FloorplanConfig(
        port_pre_assignments=port_pre_assignments,
        dse_range_min=0.6,
        dse_range_max=0.8,
    )
    floorplan_config.save_to_file(FLOOR_PLAN_CONFIG)

    left_ch = ch if ch < 8 else 8
    right_ch = 0 if ch < 8 else ch - left_ch 

    factory = get_u280_vitis_device_factory(VITIS_PLATFORM)
    #Reserving LUTs/FFs for HBM memory sub-system
    factory.reduce_slot_area(0, 0, lut=5000*2*left_ch, ff=6500*2*left_ch)
    factory.reduce_slot_area(1, 0, lut=5000*2*right_ch, ff=6500*2*right_ch)
    #Excluding DSPs on the boundary between dynamic/static region
    factory.reduce_slot_area(1, 1, dsp=100)
    factory.generate_virtual_device(Path(DEVICE_CONFIG))

    #get_u280_precollected_device(DEVICE_CONFIG)

    config = PipelineConfig(
        pp_scheme="single",
        pipeline_data_of_intra_slot_fifo="false",
        max_seconds="60",
        max_workers="1",
    )
    config.save_to_file(Path(PIPELINE_CONFIG))


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate configuration files.")
    parser.add_argument(
        "--ch",
        type=int,
        required=True,
        help="The number of in/out channels used by NTT.",
    )
    args = parser.parse_args()

    # Call the function with the provided ch
    gen_config(ch=args.ch)
