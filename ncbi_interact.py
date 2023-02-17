#!/usr/bin/env python3
"""A main class for preparing data and interacting with ncbi.
"""

from modules.ncbi import NCBI
from modules.arguments import add_arguments
import argparse
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    add_arguments(p)
    args = p.parse_args()

    # TODO: add this to NCBI init
    if args.test_dir and "test" not in args.outdir.name:
        print("changing outdir --> test")
        args.outdir:Path = args.outdir.parent / f"{args.outdir.name}_test"
    return args

def main():
    print("args:")
    args = parse_args()
    print(args)
    ncbi = NCBI(
        fastq_dir = args.fastq_dir,
        seq_report = args.seq_report,
        plate = args.plate,
        outdir = args.outdir,
        config = args.config,
        gisaid_log = args.gisaid_log,
        primer_map = args.primer_map,
        primer_scheme = args.primer_scheme,
        fasta = getattr(args,"fasta",None),
        ncbiUser = getattr(args,"username",None),
        ncbiPass = getattr(args,"password",None),
        host = getattr(args,"host",None),
        test_dir = args.test_dir,
        test_mode = getattr(args,"test_mode",False),
        )

    print("Preparing data")
    ncbi.write_presubmission_metadata()

if __name__ == "__main__":
    main()
