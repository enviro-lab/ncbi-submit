#!/usr/bin/env python3
"""A main class for preparing data and interacting with ncbi.
"""
#Git Push test. 
from ncbi_submit.ncbi import NCBI
from ncbi_submit.arguments import parse_args
from ncbi_submit.helpers import remove_empty_file,samplesFromSpuidSpecifications
import logging
logging.getLogger().addHandler(logging.StreamHandler())

def main():
    """Do NCBI submission preparation, submission, or checking, based on command line arguments"""

    args = parse_args()
    # print(args.update_xml)
    # exit(1)

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
            plate = getattr(args,"plate",None),
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
            allow_submitted = getattr(args,"update_reads",False) or getattr(args,"update_xml",False),
            update_xml = getattr(args,"update_xml",False),
            )
        spuid_endings = getattr(args,"spuid_endings",None)
        if spuid_endings: spuid_endings = spuid_endings.strip("'\"")
        
        report_files=getattr(args,"report_files",[]) or []

        if args.action == "file_prep":
            print("Preparing data")
            update_only = getattr(args,"update_only",False)
            updatedSamples = samplesFromSpuidSpecifications(spuid_endings) if update_only else []
            ncbi.write_presubmission_metadata(update_reads=args.update_reads,spuid_endings=spuid_endings,report_files=report_files,download_reports=args.download_reports,update_only=update_only,updatedSamples=updatedSamples)

            if args.prep_genbank:
                ncbi.write_genbank_submission_zip()

        elif args.action == "ftp":

            # print(args)
            # exit(1)

            if args.ftp_action == "submit":
                print(f"Submitting to {args.db}")
                ncbi.submit(db=args.db,attempt_num=args.attempt_num)

            elif args.ftp_action == "check":
                print(f"Checking {args.db}")
                ncbi.check(db=args.db,attempt_num=args.attempt_num,simple=args.simple)

            elif args.ftp_action == "get-accessions":
                print(f"Getting accessions for all reports from {args.db}")
                ncbi.get_all_accessions(db=args.db,report_files=report_files,download_reports=args.download_reports)

        elif args.action == "example":
            from example import file_getter
            file_getter.get_files(outdir=args.outdir,config=args.config,template=args.template)
            
    # remove logfile if not used
    remove_empty_file(args.logfile)

if __name__ == "__main__":
    main()
