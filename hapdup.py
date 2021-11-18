#!/usr/bin/env python3


"""
This script sets up environment paths
and invokes HapDup
"""

import os
import sys

if sys.version_info < (3,):
        raise SystemExit("Requires Python 3")

#Setting environment for local run without installation
hap_dup_root = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, hap_dup_root)

os.environ["PEPPER_MODEL_DIR"] = os.path.join(hap_dup_root, "pepper_models")
os.environ["MARGIN_CONFIG_DIR"] = os.path.join(hap_dup_root, "submodules", "margin", "params", "phase")

#entry point
from hapdup.main import main
sys.exit(main())
