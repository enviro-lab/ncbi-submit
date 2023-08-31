#!/usr/bin/env python3
from pathlib import Path
import textwrap, argparse

def add_file_prep_args(parser_file_prep:argparse.ArgumentParser):
    """Adds arguments for file preparation"""
    parser_file_prep.add_argument("--fastq_dir",type=Path,required=True,
        help="Path to directory containing one compiled fastq file for each sample")
    parser_file_prep.add_argument("--seq_report",type=Path,required=True,
        help="Path to sequencing report CSV")
    parser_file_prep.add_argument("--plate",required=True,
        help="Unique run or plate identifier")
    parser_file_prep.add_argument("--barcode_map",type=Path,required=False,default=None,
        help="Path to barcode map TSV. Used as reference to check that all samples are accounted for.")
    parser_file_prep.add_argument("--primer_map",type=Path,required=False,
        help="Path to TSV mapping of ['sample_name','Seq ID',or 'Sample #'] to ['Primer Scheme']. Used to set scheme and protocol details based on values in config.")
    parser_file_prep.add_argument("--primer_scheme",type=str,required=False,
        help="Default primer scheme for all samples. Only used if `--primer_map` not provided or contains no 'Primer Scheme' column. Used to set scheme and protocol details based on values in config.")
    parser_file_prep.add_argument("--outdir",type=Path,required=True,
        help="Path to output directory (default = ./ncbi)",default=Path("ncbi"))
    parser_file_prep.add_argument("--controls",
        help="list of names used to reference control samples (seperated by '|')")
    parser_file_prep.add_argument("--config",default=None,
        help="Path to config file"),
    # parser_file_prep.add_argument("--template",required=False,type=Path,
    #     help="Path to submission template file. If present, this overrides the template path from the config. Can be created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/")
    parser_file_prep.add_argument("--gisaid_log",required=False,type=Path,
        help="Path to submission log from gisaid upload (or other filter file in the format `virus_name; accession`). Used to indicate samples that are good enough for submission and to add gisaid accession to GenBank submission.")
    parser_file_prep.add_argument("--fasta",required=False,type=Path,
        help="Path to consensus fasta")
    parser_file_prep.add_argument("--test_dir",action="store_true",
        help="Sets submission directory to submit/Test. If not indicated, submission directory is submit/Production")
    # parser_file_prep.add_argument("--sra_only",action='store_true',
    #     help="Indicates that only the SRA tsv needs to be produced") #    # TODO: remove                 (NOTE: still testing)
    parser_file_prep.add_argument("--biosample_accessions",type=Path,required=False,
        help="Path to biosample accession file (for mapping accessions to samples)")
    parser_file_prep.add_argument("--prep_genbank",action="store_true",
        help="Indicates that GenBank XML and Zip files should be produced, adding in BioSample accessions if available.")
    parser_file_prep.add_argument("--use_existing",action="store_true",
        help="Indicates to use existing TSV files (biosample_attributes.tsv, genbank-metadata.tsv, and sra-metadata.tsv) rather than creating them all anew. Any XML or Zip submission files will be created fresh from those TSVs. If the TSVs don't exist, they will be created.")
    parser_file_prep.add_argument("--vary_spuid",action="store_true",
        help="Indicates to use append `attempt_num` to SPUID. Sometimes SPUID for updates needs to be different than SPUID from older submission.")
    parser_file_prep.add_argument("--update_reads",action="store_true",
        help="Seeks previously uploaded sample accessions. XML sheets will only include the SRA action for samples that are in the metadata but have existing submissions. Other samples will have both action sections like normal.")
    parser_file_prep.add_argument("--spuid_endings",
        help="Adds unique, explicit suffix to SPUIDs for specified samples, e.g. 'suffix1:samp1,samp2;suffix2:samp3'. Useful if updating reads from previous submissions. If a single value is given, all samples will receive that suffix.")
    parser_file_prep.add_argument("-f","--report_files",type=Path,required=False,nargs="*",default=[],
        help="Path(s) to report file(s) from which to retrieve accessions. Only used if `--update_reads` is specified. If not provided, reports will be downloaded to outdir/reports_%(bioproject_accession)")
    parser_file_prep.add_argument("-d","--download_reports",action='store_true',
        help="A flag to download report*.xml files from NCBI to outdir/reports_%(bioproject_accession). Only used if `--update_reads` is specified.")
    parser_file_prep.add_argument("--update_xml",action="store_true",
        help="Flag to only update files, not submit any reads, and allow previously submitted files to be included.")

def add_common_ftp_args(parser:argparse.ArgumentParser):
    """Adds arguments to subparser that should be the same for all ftp actions"""
    parser.add_argument("-u","--username",type=str,required=False,
        help="username - This can also be retrieved from the environmental variable 'ncbiUser'")
    parser.add_argument("-p","--password",type=str,required=False,
        help="password - This can also be retrieved from the environmental variable 'ncbiPass'")
    parser.add_argument("--host",type=str,default="ftp-private.ncbi.nlm.nih.gov",
        help="host address")
    parser.add_argument("--config",default=None,
        help="Path to config file")
    parser.add_argument("--test_dir",action="store_true",
        help="Sets submission directory to submit/Test. If not indicated, submission directory is submit/Production")

def add_submit_check_ftp_args(parser:argparse.ArgumentParser):
    """Adds arguments to subparser that should be the same for submit and check ftp actions"""
    parser.add_argument("-o","--outdir",type=Path,required=True,
        help="Path to output directory (default = ./ncbi)",default=Path("ncbi")),
    parser.add_argument("-s","--subdir",type=str,required=False,
        help="Prefix for remote subdirectory where files will be uploaded. `db` and `attempt_num` will be added to form the full subdirectory name. If not provided, defaults to `--plate`")
    parser.add_argument("--plate",required=True,
        help="Unique run or plate identifier")
    parser.add_argument("-c","--controls",type=str,
        help="'|'-delimited list of names used as controls (for excluding related files from upload)",default="Pos|Neg|NTC|PC|NEC")
    parser.add_argument("-n","--attempt_num",type=int,
        help="the attempt number of a run - used to make unique directories if multiple pre-submission directories have to be used in NCBI (errors with previous submissions)",default=1)

def add_ftp_args(parser_ftp:argparse.ArgumentParser):
    """Adds arguments for ftp interactions"""
    ftp_subparsers = parser_ftp.add_subparsers(help="Interact with ncbi via ftp",dest="ftp_action")
    ftp_subparsers.required = True

    parser_ftp_submit = ftp_subparsers.add_parser("submit",help="(for sra and biosample) - submits sra_biosample.xml as submission.xml, any referenced fastqs, and submit")
    add_common_ftp_args(parser_ftp_submit)
    # parser_ftp.add_argument("--submit",action="store_true",help="(for sra and biosample) - submits sra_biosample.xml as submission.xml, any referenced fastqs, and submit.ready")
    parser_ftp_submit.add_argument("--db",choices=["bs_sra","gb"],
        help="The ncbi submission database which determines the submisssion files as follows:"
        "   bs_sra          * sra_biosample.xml as submission.xml"
        "   (BioSample,sra) * any referenced fastqs"
        "   gb              * genbank.zip"
        "   (GenBank)       * genbank.xml as submission.xml"
        )
    parser_ftp_submit.add_argument("-f","--fastq_dir",type=Path,required=True,
        help="Path to local directory containing one compiled fastq file for each sample")
    parser_ftp_submit.add_argument("-T","--test_mode",action="store_true",
        help="Everything but the upload will happen (including signing into and navigating the ftp site). Files that would be transfered will be listed.")
    parser_ftp_submit.add_argument("--update_xml",action="store_true",
        help="Flag to only update files, not submit any reads, and allow previously submitted files to be included.")
    add_submit_check_ftp_args(parser_ftp_submit)
    # parser_ftp_submit.add_argument("--sra_only",action='store_true',
    #     help="Indicates that only the SRA tsv needs to be produced") #                                      (NOTE: still testing)

    parser_ftp_check = ftp_subparsers.add_parser("check",help="Reports on submission success. If no `db` specified, all will be checked")
    add_common_ftp_args(parser_ftp_check)
    parser_ftp_check.add_argument("--db",choices=["bs_sra","gb"],
        help="The ncbi database for which to check submission status")
    # parser_ftp.add_argument("--check",action="store_true",
    #     help="Report on submission success. If no `db` specified, all will be checked")
    parser_ftp_check.add_argument("--simple",action="store_true",
        help="Limits submission check to specified submission `db` and `attempt_num`. If not specified, all possible `attemp_num`s will be checked")
    add_submit_check_ftp_args(parser_ftp_check)
    # parser_ftp_check.add_argument("-o","--outdir",type=Path,required=True,
    #     help="Path to output directory (default = ./ncbi)",default=Path("ncbi")),
    # parser_ftp_check.add_argument("-s","--subdir",type=str,required=False,
    #     help="Prefix for remote subdirectory where files will be uploaded. `db` and `attempt_num` will be added to form the full subdirectory name. If not provided, defaults to `--plate`")
    # parser_ftp_check.add_argument("--plate",required=True,
    #     help="Unique run or plate identifier")
    # parser_ftp_check.add_argument("-c","--controls",type=str,
    #     help="'|'-delimited list of names used as controls (for excluding related files from upload)",default="Pos|Neg|NTC|PC|NEC")
    # parser_ftp_check.add_argument("-n","--attempt_num",type=int,
    #     help="the attempt number of a run - used to make unique directories if multiple pre-submission directories have to be used in NCBI (errors with previous submissions)",default=1)

    parser_ftp_get_accessions = ftp_subparsers.add_parser("get-accessions",help="Creates accession map from all current report files in NCBI for `db`")
    add_common_ftp_args(parser_ftp_get_accessions)
    parser_ftp_get_accessions.add_argument("--db",choices=["bs_sra","gb"],
        help="The ncbi database for which to download report files and retrieve accessions")
    parser_ftp_get_accessions.add_argument("-o","--outdir",type=Path,required=False,default=None,
        help="Path to output directory (default: `reports_dir` from config)")
    parser_ftp_get_accessions.add_argument("-f","--report_files",type=Path,required=False,nargs="*",
        help="Path(s) to report file(s) from which to retrieve accessions. If not provided, reports will be downloaded to outdir/reports_%(bioproject_accession)")
    parser_ftp_get_accessions.add_argument("-d","--download_reports",action='store_true',
        help="A flag to download report*.xml files from NCBI to outdir/reports_%(bioproject_accession)")
    # parser_ftp.add_argument("--get-accessions",action="store_true",
    #     help="Creates accession map from all current report files in NCBI for `db`")
    

    # parser_ftp.add_argument("--biosample_accessions",type=Path,required=False,
    #     help="Path to existing biosample accession TSV (mapping accessions to samples)")
    # parser_ftp.add_argument("--gisaid_log",required=True,type=Path,
    #     help="Path to submission log from gisaid upload - used to weed out unwanted seqs from `--fasta`")

def add_example_args(parser_example:argparse.ArgumentParser):
    """Adds arguments for acquiring examples"""
    parser_example.add_argument("--config",action="store_true",help="Only write the example config file")
    parser_example.add_argument("--template",action="store_true",help="Only write the example template file")
    parser_example.add_argument("--outdir",required=False,type=Path,help="Path to output directory (default = ./ncbi)",default=Path("ncbi"))

def add_arguments(parser:argparse.ArgumentParser):
    """Add all aruments to ncbi_interact argument parser

    Args:
        `parser`: ArgumentParser
    """

    from ncbi_submit.version import __version__
    # parser.add_argument("--plate",required=True,help="Unique run or plate identifier")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s ("+__version__+")")
    parser.add_argument("--log",choices=["DEBUG","WARNING","INFO"],required=False,default="",help="Sets logging level.")
    subparsers = parser.add_subparsers(
        help="Prepare NCBI submission TSVs or submit to NCBI using one of these actions",dest="action")

    parser_file_prep = subparsers.add_parser("file_prep",formatter_class=argparse.RawTextHelpFormatter,
        help=textwrap.dedent("""\
        Create NCBI-style submission TSVs by combining 
          sequencing report metadata, primer scheme info, and EpiCov accessions. 
        Outputs include: 
            'biosample_attributes.tsv'
            'sra-metadata.tsv'
            'genbank-metadata.tsv'
            'merged-metadata.tsv' - combination of fields in above files
        """))
    
    parser_ftp = subparsers.add_parser("ftp",formatter_class=argparse.RawTextHelpFormatter,
        help=textwrap.dedent(
        """\
        NOTE: Only use one or the other of `submit` and `check`.

        Submissions/checks differ based on the flag `--db`:
        * "bs": Submits BioSample metadata.
        * "sra": Submits fastqs and SRA metadata, linking BioSample accessions if available.
        * "bs_sra": Submits BioSample and SRA metadata linked, plus fastqs. (recommended over 'bs' or 'sra' individually)
        * "gb": Submits fastas and genbank metadata, linking BioSample accessions if available.

        `submit`: Submits files to desired `--db`.
        `check`: Downloads report.xml and reports on submission status. Checks all databases if `--db` not specified.
        `get-accessions`: Creates accession map from all current report files in NCBI for `db`.
        """))
    
    parser_example = subparsers.add_parser("example",formatter_class=argparse.RawTextHelpFormatter,
        help=textwrap.dedent(
        """\
        Get templates for "config.py" and/or "template.sbt" to a specified `outdir`
        """))

    add_file_prep_args(parser_file_prep)
    add_ftp_args(parser_ftp)
    add_example_args(parser_example)

    return parser