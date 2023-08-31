#!/usr/bin/env python3
"""A set up useful tools/functions that aren't specific to the other classes
"""

from functools import partial
import gzip
import logging
from typing import Literal
import pandas as pd
from pathlib import Path
import io
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from ncbi_submit.report import Report

main_dir = list(Path(__file__).resolve().parents)[1]

class DataHolder():
    """Store filenames with associated df"""

    def __init__(self,file):
        self.file = file
        self.df = pd.DataFrame()
    def set_df(self,df):
        self.df = df

def warn(message):
    """Prints a message and exits the script with a failing return code"""
    print(message)
    exit(1)

def getConfig(config_file=None):
    """Returns `config_file` as a dictionary and it's path"""

    # decide config file (use example if not given)
    if config_file==None:
        config_file = main_dir / "config/config.py"
    else: config_file = Path(config_file)
    # ensure config exists
    if not config_file.exists():
        warn(f"{config_file} does not exist. Using default values. "\
            "`--ncbi_config` file must be provided")
    # add config vars to `config` dict and return
    config = {}
    logging.info("reading config:")
    exec(io.open(config_file).read(), config)
    return {k:v for k,v in config.items() if k!='__builtins__'}, config_file

def hasUnderscoredHeaders(columns):
    """Returns True if any field names in `columns` start with '_'"""

    return "_sample_name" in columns or "_geo_loc_name" in columns

def isBioProjectAccession(accession):
    """Returns True if `accession` looks like a BioProject accession"""

    if accession==None:
        return False
    return True in [accession.startswith(x) for x in ("PRJE","PRJN","PRJD","PSUB")]

def checkBioProjectAccessions(accessions):
    """Passes if all `accessions` look like accessions, else raises ValueError

    
        Args:
            `accession` (Collection): accessions to check

        Raises:
            ValueError: if accessions appear invalid
    """

    for acc in accessions:
        if not isBioProjectAccession(acc):
            raise ValueError(f"The provided BioProject accession '{acc}' does not fit ther standard format (PRJ[D|E|N]xxxxxx or PSUBxxxxxx).")

def finalize_df(df,final_cols=[]):
    """Returns `df` with `final_cols` as the only columns or else all columns
    
    Args:
        `df` (DataFrame): the base dataset to write out
        `final_cols` (list): if provided, only these columns will be included
    """

    return  df[final_cols] if final_cols else df.copy()

def df_2_tsv(df,outfile,name=None):
    """Writes out df as tsv file (without index)
    
    Args:
        `df` (DataFrame): the base dataset to write out
        `outfile` (str | Path): where to write TSV
        `name` (str): if provided, `name` and `outfile` are logged in stdout
    """

    if name:
        print(f"Writing `{name}` outfile:\n   {outfile}")
    outdir = Path(outfile).parent
    outdir.mkdir(parents=True,exist_ok=True)
    df.to_csv(outfile, index=False, sep="\t")

def series2dict(series:pd.Series,limit_to:list=None):
    """Converts series to list

    Args:
        `series` (pd.Series): The series to convert to a dict
        `limit_to` (list, optional): If provided, only the listed attributes will apear in the final dict
    """

    if limit_to == None:
        return series.to_dict()
    else:
        return {k:v for k,v in series.to_dict().items() if k in limit_to}

def check_missing(actual_cols,required_cols,name):
    """Verifies no required columns are missing from given database df"""

    missing = set(required_cols) - set(actual_cols)
    if len(missing) > 0:
        s = "" if len(missing) == 1 else "s"
        raise Exception(f"Required {name} column{s} not found:\n",missing)

def remove_item(unwanted_item,a_list):
    """Returns list without specified column"""

    return [col for col in a_list if col != unwanted_item]

def remove_items(unwanted_items,a_list):
    """Returns a list with each item in `unwanted_items` excluded from `a_list`"""

    return list(set(a_list)-set(unwanted_items))

def ensureAllSamplesHaveNames(biosample:pd.DataFrame):
    """If any sample_name values are NaN in `biosample` DataFrame, print this warning and exit"""

    missing:pd.DataFrame = biosample[biosample["sample_name"].isna()]
    if missing.empty: return
    print('\n\n\nWARNING:\n')
    print("The following samples are missing their `sample_name`. \n"
    "This very likely means they were submitted to gisaid in an old run and added to \n"
    " the most recent gisaid submission logfile despite their barcodes not making it \n"
    " through the most recent artic pipeline run. \n")
    print(missing[['sample_name','gisaid_virus_name']])
    print("\nSteps to fix:\n"
    "If the missing samples are in a Sequencing-report-<plate>.csv file from an old run:\n"
    "   * add them to the new report\n"
    "   * run `file_prep` again without submitting\n"
    "   * verify this message doesn't show up again and that missing samples show up\n"
    "   * if so, you should be fine to submit\n"
    "If not, check through gisaid log and other CSVs to see where things are missing... good luck\n")
    exit(1)

def get_bioproject_spuid(ncbi):
    """Returns SPUID to use if creating new BioProject
    
    Args:
        `ncbi` (NCBI): An object storing information about the submission data, goals, and useful methods
    """

    ncbi:NCBI
    return f'<SPUID spuid_namespace="{ncbi.centerAbbr}">{ncbi.bioproject_presets["spuid"]}</SPUID>'

def asDate(some_date):
    """Converts ftp.mlsd() `modify` time to easy-to-read date/time"""
    
    y,mo,d,h,m,s = some_date[:4],some_date[4:6],some_date[6:8],some_date[8:10],some_date[10:12],some_date[12:14]
    return f"{y}-{mo}-{d} {h}:{m}:{s}"

def remove_empty_file(file):
    """Removes `file` if empty"""

    file = Path(file)
    if not file.exists(): return
    lines = file.open('r').readlines()
    if len(lines) == 0:
        file.unlink()
        return
    if len(lines) > 1: return
    else:
        line = lines[0]
        if len(line) > 1:
            return
        if line == "\n" or not line:
            file.unlink()

def ensure_outdir_viable(outdir:Path):
    """Ensures `outdir` isn't a file and makes the directory, if needed
    
    Returns:
        Path(outdir)
    """

    if not type(outdir) == Path:
        outdir = Path(outdir).resolve()
    if outdir.is_dir(): return outdir
    if not outdir.exists():
        outdir.mkdir(exist_ok=True,parents=False)
        return outdir
    if outdir.is_file(): raise FileExistsError(f"Proposed `outdir` ({outdir}) already exists as a file.")

def is_fastq(file):
    """Returns True if file is fastq or fastq.gz else False

    Args:
        file (str|Path): filename to check
    """

def bioproject_in_report(report,bioproject):
    """Returns True if `bioproject` is in `report`"""

    with open(report) as fh:
        for line in fh:
            if bioproject in line:
                return True
    return False

def get_accession_dict_from_file(report_file,db:Literal["bs_sra","bs","sra","gb",None]=None,bioproject=None,ignore_failure=True):
    """Returns dict of accessions for each sample in `report_file`
    
    Args:
        `report_file` (str | Path): path to NCBI report file to read in
        `db` (str): name of database for which to get accessions
        `bioproject` (str): BioProject accession to include exclusively, if provided
    """

    if bioproject != None and not bioproject_in_report(report_file,bioproject):
        return {}
    report = Report(report_file)
    accession_dict = {}
    if db == "bs_sra":
        biosample_accessions = report.getAccessionDict("BioSample",by_sample_name=True,ignore_failure=ignore_failure)
        sra_accessions = report.getAccessionDict("SRA",by_sample_name=True,ignore_failure=ignore_failure)
        for sample_id in set(biosample_accessions.keys()).union(sra_accessions.keys()):
            biosample_accession = biosample_accessions.get(sample_id)
            sra_accession = sra_accessions.get(sample_id)
            accession_dict[sample_id] = {"BioSample":biosample_accession,"SRA":sra_accession}
    elif db == "bs":
        biosample_accessions = report.getAccessionDict("BioSample",by_sample_name=True,ignore_failure=ignore_failure)
        for sample_id,biosample_accession in biosample_accessions.items():
            accession_dict[sample_id] = {"BioSample":biosample_accession}
    elif db == "sra":
        sra_accessions = report.getAccessionDict("SRA",by_sample_name=True,ignore_failure=ignore_failure)
        for sample_id,sra_accession in sra_accessions.items():
            accession_dict[sample_id] = {"SRA":sra_accession}
    else:
        raise NotImplementedError("Not yet prepared to handle databases other than 'bs_sra','bs','sra'")
    return accession_dict

def is_fastq(file):
    """Returns True if file is fastq or fastq.gz else False

    Args:
        file (str|Path): filename to check
    """

    file = Path(file)
    if "fastq" not in file.suffixes:
        return False
    _open = partial(gzip.open,mode="rt") if file.suffix == "gz" else open
    with _open(file) as fh:
        for title, seq, qual in FastqGeneralIterator(fh):
            if title and seq and qual:
                return True
            else: return False