#!/usr/bin/env python3
"""A main class for preparing data and interacting with ncbi.
"""
#Git Push test. 
from ncbi_submit.ncbi import NCBI
from ncbi_submit.arguments import add_arguments
from ncbi_submit.helpers import remove_empty_file,ensure_outdir_viable
import argparse, logging

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    add_arguments(p)
    args = p.parse_args()
    # set logging level, if needed
    setattr(args,"logfile",args.outdir/'ncbi-submit.log')
    args.logfile.parent.mkdir(exist_ok=True,parents=True)
    args.outdir = ensure_outdir_viable(args.outdir)
    logging.basicConfig(filename=args.logfile, encoding='utf-8', level=getattr(logging,args.log,None))
    return args

def main():
    args = parse_args()

    if args.action == "example":
        from example import file_getter
        # if neither option is passed, assume both files are desired
        if not args.config and not args.template: args.config,args.template=True,True
        # write files
        file_getter.get_files(outdir=args.outdir,config=args.config,template=args.template)
    
    else:

        logging.info("args:")
        logging.info(args)

        ncbi = NCBI(
            fastq_dir = getattr(args,"fastq_dir",None),
            seq_report = getattr(args,"seq_report",None),
            barcode_map = getattr(args,"barcode_map",None),
            plate = args.plate,
            subdir = getattr(args,"subdir",None),
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

        elif args.action == "ftp":

            if args.submit:
                print(f"Submitting to {args.db}")
                ncbi.submit(db=args.db,attempt_num=args.attempt_num)

            if args.check:
                print(f"Checking {args.db}")
                ncbi.check(db=args.db,attempt_num=args.attempt_num,simple=args.simple)

        elif args.action == "example":
            from example import examples
            file_getter.get_files(outdir=args.outdir,config=args.config,template=args.template)
            
    # remove logfile if not used
    remove_empty_file(args.logfile)

if __name__ == "__main__":
    main()
