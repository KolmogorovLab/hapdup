#!/usr/bin/env python3


"""
This script sets up environment paths
and invokes viralFlye without installation.
"""

import os
import sys

if sys.version_info < (3,):
        raise SystemExit("Requires Python 3")

#Setting environment for local run without installation
hap_dup_root = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, hap_dup_root)

#entry point
from hap_dup.main import main
sys.exit(main())
