"""Generate rapidstream configuration files."""
__copyright__ = """
Copyright (c) 2024 RapidStream Design Automation, Inc. and contributors.
All rights reserved. The contributor(s) of this file has/have agreed to the
RapidStream Contributor License Agreement.
"""
import os
from pathlib import Path
import argparse

from rapidstream.assets.device_library.u280.u280 import get_u280_precollected_device
from rapidstream.assets.floorplan.floorplan_config import FloorplanConfig

TEMP_DIR = Path("config_build")
FLOOR_PLAN_CONFIG = TEMP_DIR / "floorplan_config.json"
DEVICE_CONFIG = TEMP_DIR / "device_config.json"


def gen_config(num_ch: int) -> None:
    """Generate configuration files."""
    if not os.path.exists(TEMP_DIR):
        os.makedirs(TEMP_DIR)

    port_pre_assignments = {
        "ap_clk": "SLOT_X0Y0:SLOT_X1Y2",
        "ap_rst_n": "SLOT_X0Y0:SLOT_X1Y2",
        "interrupt": "SLOT_X0Y0:SLOT_X1Y2",
        "s_axi_control_.*": "SLOT_X0Y0:SLOT_X1Y2",
    }
    for i in range(num_ch):
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
    get_u280_precollected_device(DEVICE_CONFIG)


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate configuration files.")
    parser.add_argument(
        "--num_ch",
        type=int,
        required=True,
        help="The number of in/out channels used by NTT.",
    )
    args = parser.parse_args()

    # Call the function with the provided num_ch
    gen_config(num_ch=args.num_ch)
