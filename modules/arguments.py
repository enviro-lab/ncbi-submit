#!/usr/bin/env python3
from pathlib import Path
import textwrap, argparse

def add_file_prep_args(parser_file_prep):
    ## file_prep
    parser_file_prep.add_argument("--fastq_dir",type=Path,required=True,
        help="Path to directory containing one compiled fastq file for each sample")
    parser_file_prep.add_argument("--seq_report",type=Path,required=True,
        help="Path to sequencing report CSV")
    parser_file_prep.add_argument("--plate",required=True,
        help="Unique run or plate identifier")
    parser_file_prep.add_argument("--primer_map",type=Path,required=False,
        help="Path to TSV mapping of ['sample_name','Seq ID',or 'Sample #'] to ['Primer Scheme']. Used to set scheme and protocol details based on values in config.")
    parser_file_prep.add_argument("--primer_scheme",type=str,required=False,
        help="Only used if `--primer_map` not provided. Used to set scheme and protocol details based on values in config.")
    parser_file_prep.add_argument("--outdir",type=Path,required=True,
        help="Path to output directory (default = ./ncbi)",default=Path("ncbi"))
    parser_file_prep.add_argument("--controls",
        help="list of names used to reference control samples (seperated by '|')")
    # parser_file_prep.add_argument("--protocol",required=True,
        # help="'sequencing_protocol_name' for this batch of samples")
    # parser_file_prep.add_argument("--email",type=str,required=False, # TODO: remove
    #     help="sequence submitter email address")
    # parser_file_prep.add_argument("--center",type=str,required=True,dest='center_abbr',
    #     help="center abbreviation for NCBI programmatic submissions")
    parser_file_prep.add_argument("--config",default=None,
        help="Path to multisub config file"),
    parser_file_prep.add_argument("--submitted",default=None,
        help="Path to multisub config file"),
    # parser_file_prep.add_argument("--template",required=False,type=Path,
    #     help="Path to submission template file. If present, this overrides the template path from the config. Can be created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/")
    parser_file_prep.add_argument("--gisaid_log",required=False,type=Path,
        help="Path to submission log from gisaid upload (of filter file in the format `virus_name; accession`) - used to indicate allowed seqs from `--fasta` and to add gisaid accession to GenBank submission")
    # parser_file_prep.add_argument("--fasta",required=True,type=Path,
    #     help="Path to consensus fasta")
    parser_file_prep.add_argument("--test_dir",action="store_true",
        help="Sets submission directory to submit/Test. If not indicated, submission directory is submit/Production")
    parser_file_prep.add_argument("--sra_only",action='store_true',
        help="Indicates that only the SRA tsv needs to be produced") #                                      (NOTE: still testing)
    parser_file_prep.add_argument("--biosample_accessions",type=Path,required=False,
        help="Path to biosample accession file (for mapping accessions to samples)")

def add_ftp_args(parser_ftp):
    ## ftp
    parser_ftp.add_argument("--db",help="options: [bs_sra, gb]"
        "  if `--submit` is selected the following files are submitted"
        "     bs_sra          * sra_biosample.xml as submission.xml"
        "     (BioSample,sra) * any referenced fastqs"
        "     gb              * genbank.zip"
        "     (GenBank)       * genbank.xml as submission.xml"
        "  if `--check` is selected, the corresponding report files are downloaded and checked")
    parser_ftp.add_argument("--submit",action="store_true",help="(for sra and biosample) - submits sra_biosample.xml as submission.xml, any referenced fastqs, and submit.ready")
    # parser_ftp.add_argument("--submit_bs_sra",action="store_true",help="(for sra and biosample) - submits sra_biosample.xml as submission.xml, any referenced fastqs, and submit.ready")
    # parser_ftp.add_argument("--submit_gb",action="store_true",help="(for genbank) - submits genbank.zip, genbank.xml as submission.xml, and submit.ready")
    # parser_ftp.add_argument("--check",help="options: [bs_sra, gb]. Looks for report file related to given submission. Downloads to outdir/bs_sra_reports or outdir/gb_reports")
    parser_ftp.add_argument("--check",action="store_true",help="verify successful submission or report on issues")
    parser_ftp.add_argument("--simple",action="store_true",help="if checking on submission, this option limits output to one line indicating submission status")
    parser_ftp.add_argument("-u","--username",type=str,required=True,help="username")
    parser_ftp.add_argument("-p","--password",type=str,required=True,help="password")
    parser_ftp.add_argument("-o","--host",type=str,required=True,help="host address")
    parser_ftp.add_argument("--outdir",type=Path,required=True,
        help="Path to output directory (default = ./ncbi)",default=Path("ncbi")),
    parser_ftp.add_argument("-s","--subdir",type=str,required=True,
        help="remote subdirectory where files will be uploaded")
    parser_ftp.add_argument("-f","--fastq_dir",type=Path,required=True,
        help="Path to local directory containing one compiled fastq file for each sample")
    parser_ftp.add_argument("-c","--controls",type=str,
        help="'|'-delimited list of names used as controls (for excluding related files from upload)",default="Pos|Neg|NTC|PC|NEC")
    parser_ftp.add_argument("-n","--attempt_num",type=int,
        help="the attempt number of a run - used to make unique directories if multiple pre-submission directories have to be used in NCBI (errors with previous submissions)",default=1)
    # parser_ftp.add_argument("-a","--ask_about_possible_controls",action="store_true",
    #     help="If selected, you'll be asked whether files assumed to be controls should still be uploaded. (Default=False)",default=False)
    parser_ftp.add_argument("-T","--test_mode",action="store_true",
        help="Everything but the upload will happen (including signing into and navigating the ftp site). Files that would be transfered will be listed.")
    parser_ftp.add_argument("--config",default=None,
        help="Path to multisub config file")
    parser_ftp.add_argument("--test_dir",action="store_true",
        help="Sets submission directory to submit/Test. If not indicated, submission directory is submit/Production")
    parser_ftp.add_argument("--sra_only",action='store_true',
        help="Indicates that only the SRA tsv needs to be produced") #                                      (NOTE: still testing)
    parser_ftp.add_argument("--biosample_accessions",type=Path,required=False,
        help="Path to biosample accession file (for mapping accessions to samples)")
    # parser_ftp.add_argument("--gisaid_log",required=True,type=Path,
    #     help="Path to submission log from gisaid upload - used to weed out unwanted seqs from `--fasta`")

def add_biosample_args(parser_biosample):
    ## add_biosample
    parser_biosample.add_argument("--outdir",type=Path,default=Path("ncbi"),required=True,
        help="Path to output directory (default = ./ncbi). Must contain the file 'genbank-metadata.tsv'")
    parser_biosample.add_argument("--config",default=None,
        help="Path to multisub config file")
    parser_biosample.add_argument("--fasta",required=True,type=Path,
        help="Path to consensus fasta containing all desired samples. Extras are allowed, as this will be filtered by either `--gisaid_log` or the 'GenBank.tsv'")
    parser_biosample.add_argument("--plate",required=True,
        help="Unique run or plate identifier")
    parser_biosample.add_argument("--sra_only",action='store_true',
        help="Indicates that only the SRA tsv needs to be produced") #                                      (NOTE: still testing)
    parser_biosample.add_argument("--biosample_accessions",type=Path,required=False,
        help="Path to biosample accession file (for mapping accessions to samples)")
    parser_biosample.add_argument("--test_dir",action="store_true",required=False,
        help="Sets submission directory to submit/Test. If not indicated, submission directory is submit/Production")
    parser_biosample.add_argument("--gisaid_log",required=False,type=Path,default=None,
        help="Path to submission log from gisaid upload - used to weed out unwanted seqs from `--fasta`")
    parser_biosample.add_argument("--template",required=False,type=Path,
        help="Path to submission template file. If present, this overrides the template path from the config. Can be created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/")


def add_arguments(parser):
    """Add all aruments to ncbi_interact argument parser

    Args:
        parser: ArgumentParser
    """

    parser.add_argument("--plate",nargs=1,default="",help="unique run or plate identifier")
    subparsers = parser.add_subparsers(
        help="Prepare NCBI submission TSVs or submit to NCBI using one of these actions",dest="action")
    parser_file_prep = subparsers.add_parser("file_prep",formatter_class=argparse.RawTextHelpFormatter,
        help=textwrap.dedent(
        "Create NCBI-style submission TSVs by combining \
        \n  sequencing report metadata, primer scheme info, and EpiCov accessions. \
        \nOutputs include: \
        \n    'biosample_attributes.tsv'\
        \n    'sra-metadata.tsv'\
        \n    'genbank-metadata.tsv'\
        \n    'merged-metadata.tsv' - combination of fields in above files\
        \n"))
    
    parser_ftp = subparsers.add_parser("ftp",formatter_class=argparse.RawTextHelpFormatter,
        help=textwrap.dedent(
        """\
        `--submit` will submit fastqs for all samples with gisaid accessions.
        Fastqs will be uploaded to SRA unless the --test_mode flag is set
        or
        `--check` will download report.xml and look for issues
        NOTE: only one or the other of `--submit` and `--check` can be used."""))
    
    parser_biosample =  subparsers.add_parser("add_biosamples",
        help="Add BioSample accessions to the previously created 'genbank-metadata.tsv'")
    
    add_file_prep_args(parser_file_prep)
    add_ftp_args(parser_ftp)
    add_biosample_args(parser_biosample)

    return parser