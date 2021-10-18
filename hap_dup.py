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

os.environ["PEPPER_MODEL"] = os.path.join(hap_dup_root, "pepper_models", "PEPPER_VARIANT_R941_ONT_V5.pkl")
os.environ["MARGIN_MODEL"] = os.path.join(hap_dup_root, "submodules", "margin", "params", "misc",
                                          "allParams.ont_haplotag.sv.json")

#entry point
from hap_dup.main import main
sys.exit(main())
