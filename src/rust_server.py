from pycpm2 import sorting
from CPMMD.parameter_file import load_file
import sys
from dataclasses import asdict


parameters = load_file(sys.argv[1]).cpm


sorting(
        parameters
    )