#!/usr/bin/env python3
import os
import sys

from acpype.cli import init_main
from acpype.params import binaries


def main():
    antechamber_bin = os.environ.get("ANTECHAMBER_BIN", "").strip()
    if antechamber_bin:
        binaries["ac_bin"] = antechamber_bin
    init_main(argv=sys.argv[1:])


if __name__ == "__main__":
    main()
