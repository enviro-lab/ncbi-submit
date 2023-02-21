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
    return args

def main():
    args = parse_args()
    print("args:")
    print(args)
    ncbi = NCBI(
        fastq_dir = getattr(args,"fastq_dir",None),
        seq_report = getattr(args,"seq_report",None),
        plate = args.plate,
        outdir = args.outdir,
        config = args.config,
        gisaid_log = getattr(args,"gisaid_log",None),
        primer_map =  getattr(args,"primer_map",None),
        primer_scheme = getattr(args,"primer_scheme",None),
        fasta = getattr(args,"fasta",None),
        ncbiUser = getattr(args,"username",None),
        ncbiPass = getattr(args,"password",None),
        host = getattr(args,"host",None),
        test_dir = args.test_dir,
        test_mode = getattr(args,"test_mode",False),
        use_existing = getattr(args,"use_existing",None),
        )

    if args.action == "file_prep":
        print("Preparing data")

        ncbi.write_presubmission_metadata()

        if args.prep_genbank:
            ncbi.write_genbank_submission_zip()


if __name__ == "__main__":
    main()
