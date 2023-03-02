#!/usr/bin/env python3
"""A main class for preparing data and interacting with ncbi.
"""

from ftplib import FTP
from io import StringIO
import os, logging
from textwrap import dedent
import pandas as pd
from ncbi_submit.helpers import *
from ncbi_submit.report import Report
from ncbi_submit.xml_format import SRA_BioSample_Submission,GenBank_Submission
from zipfile import ZipFile
from Bio import SeqIO

### TRY to use these as parameters from ncbi_interact (avoid them as attributes)
# self.sra_only = sra_only
# self.biosample_accessions = biosample_accessions
# self.subdir = subdir
# self.db = db
# self.test_mode = test_mode
# self.host = host
# self.username = username
# self.password = password
# self.submit = submit
# self.check = check
# self.attempt_num = attempt_num

class NCBI:
    def __init__(self,plate,fastq_dir,seq_report,outdir,barcode_map=None,
                 config=None,controls=None,gisaid_log=None,primer_map=None,
                 primer_scheme=None,fasta=None,ncbiUser=None,ncbiPass=None,
                 host=None,test_dir=False,test_mode=False,use_existing=False,
                 subdir=None,vary_spuid=False) -> None:
        """Creates an object that can prepare, submit, and track csv/xml submission files 

        Args:
            config: (str, Path, None) Location of the config file to use
        """

        # get any variables from config
        cf,self.config_file = getConfig(config)
        self.template = cf["template"]
        self.host = host or cf["host"] or "ftp-private.ncbi.nlm.nih.gov"
        self.ncbiUser = ncbiUser or cf.get("ncbiUser")
        self.ncbiPass = ncbiPass or cf.get("ncbiPass")
        self._check_login()
        self.centerAbbr = cf["centerAbbr"]
        self.contact = cf["contact"]
        self.email = cf["email"]
        self.phone = cf["phone"] # not used (but could be for template generation)
        self.controls = controls or cf.get("controls")
        self.ncbiCountry = cf["ncbiCountry"]
        self.affiliation = cf["affiliation"]
        self.authors = cf["authors"]
        self.bioproject_presets = cf["bioproject"]
        self.biosample_presets = cf["biosample"]
        self.genbank_presets = cf.get("genbank",{"comment":{"data":{},"include_comment":False}})
        self.sars_cov_2_diag_pcr_ct_values = cf.get("sars_cov_2_diag_pcr_ct_values",{})
        self.sra_presets = cf["sra"]
        self.allowed_schemes = cf["allowed_schemes"]
        self.protocol_scheme = cf["protocol_scheme"]
        self.scheme = cf["scheme"]
        self.submitted_samples_file = cf.get("submitted_samples")
        self.ftp = None
        self.test_dir = test_dir if test_dir==True else cf.get("test_dir",test_dir)
        self.test_mode = test_mode if test_mode==True else cf.get("test_mode",test_mode)
        self.use_existing = use_existing
        self.vary_spuid = vary_spuid

        # set plate-specific details
        self.fastq_dir = Path(fastq_dir) if fastq_dir else None
        self.seq_report = Path(seq_report) if seq_report else None
        self.plate = plate
        self.outdir = self._set_outdir(outdir)
        self.exclude_file = self.outdir / "samples2exclude.txt"
        self.barcode_map = None if barcode_map==None else Path(barcode_map) # only required for file prep
        self.gisaid_log = gisaid_log
        self.primer_map = primer_map
        self.primer_scheme = primer_scheme
        self.fasta = fasta

        # vars for ftp
        possible_report_dirs = list(self.outdir.glob("*bs*_reports")) # there should only be one (most likely)
        self.report_dir = possible_report_dirs[0] if possible_report_dirs else ""
        # self.report_dir = self.outdir / "None_reports" # TODO: remove this - used only for testing flub
        self.valid_dbs = ["sra","gb","bs_sra","bs"]
        self.subdir = self._set_default_subdir(subdir)
        self.submit_area,self.submit_dir,self.initial_submission_dirs,self.ftp_file_sizes = self._set_ftp_vars()
        self.uploaded_something = False

        self.genbank = dict(
            name = "genbank",
            tsv = self.outdir / "genbank-metadata.tsv",
            xml_text = None,
            xml_file = self.outdir / "genbank.xml",
            zip_file = self.outdir / "genbank.zip",
            df = None,
        )
        self.biosample = dict(
            name = "biosample",
            tsv = self.outdir / "biosample_attributes.tsv",
            xml_text = None,
            xml_file = self.outdir / "sra_biosample.xml",
            df = None,
        )
        self.sra = dict(
            name = "sra",
            tsv = self.outdir / "sra-metadata.tsv",
            xml_text = None,
            xml_file = self.outdir / "sra_biosample.xml",
            df = None,
        )
        self.merged = dict(
            tsv = self.outdir / "merged-metadata.tsv",
            df = None,
        )

        # dataframes that pretty much just read in input spreadsheets, if provided
        self.seq_report_df = self._get_seq_report_df()
        self.barcode_df = self._prep_barcode_df()
        self.gisaid_df = self._prep_gisaid_df()
        self.gisaid_submitted_samples = None

    ### --- vvv --- miscellaneous --- vvv --- ###

    def _set_outdir(self,outdir):
        """Insures "test" is in `outdir` basename, if using `test_dir` flag"""

        outdir = Path(outdir)
        if self.test_dir and "test" not in outdir.name:
            outdir = outdir.parent / f"{outdir.name}_test"
            print(f"Resetting outdir --> {outdir}")
            return outdir
        return outdir

    def _check_login(self,missing_ok=True):
        """Checks for login info in config or env variables

        Args:
            missing_ok (bool, optional): If False, warns about missing login details. Defaults to True.
        """

        if not self.ncbiUser:
            self.ncbiUser = os.getenv("ncbiUser")
        if not self.ncbiPass:
            self.ncbiPass = os.getenv("ncbiPass")
        if not missing_ok:
            if not self.ncbiUser or not self.ncbiPass:
                warn("`ncbiUser` and/or `ncbiPass` not provided but required for NCBI ftp interaction. They may be placed in the `ncbi_config` file or environment variables")

    def virus_to_sample_name(self,virus_name):
        """Extracts our typical sample identifier from GISAID Virus name"""

        return virus_name.split("/")[2].replace(f"{self.affiliation['sub']}-","")
    
    def get_bioproject_accession(self,default=None):
        """Returns default BioProject accession or else decides based on `test_dir` flag"""

        if isBioProjectAccession(default):
            return default
        acc = "bioproject_test_accession" if self.test_dir else "bioproject_accession"
        return self.bioproject_presets[acc]

    ### --- ^^^ --- miscellaneous --- ^^^ --- ###
    ### --- vvv --- file prep --- vvv --- ###

    def _get_seq_report_df(self):
        """Creates df from `seq_report`"""

        # local sequencing report data
        if self.seq_report == None:
            return pd.DataFrame()
            # raise Exception("`seq_report` must be provided")
        seq_report_df = pd.read_csv(self.seq_report).rename(columns={"primer_scheme":"Primer Scheme","primer scheme":"Primer Scheme",})
        # if using biosample with headers beginning with underscores
        if hasUnderscoredHeaders(seq_report_df.columns):
            new_cols = {col:col.lstrip("_") for col in seq_report_df.columns if col.startswith("_")}
            seq_report_df = seq_report_df.rename(columns=new_cols)
            seq_report_df = seq_report_df.dropna(how='all')
        else:
            # if using initial EML rather than sequencing report:
            if "Seq ID" in seq_report_df.columns:
                seq_report_df = seq_report_df[[col for col in seq_report_df.columns if col != "Sample #"]]
                seq_report_df = seq_report_df.rename(columns={"Seq ID":"Sample #"})
            # weed out controls
            if "Sample #" in seq_report_df.columns:
                seq_report_df = seq_report_df[seq_report_df["Sample #"].str.contains(self.controls)==False]
            if "Test date" in seq_report_df.columns:
                seq_report_df["Test date"] = pd.DatetimeIndex(seq_report_df['Test date']).strftime('%Y-%m-%d')
            if 'Sample ID' in seq_report_df.columns:
                seq_report_df['Sample ID'] = seq_report_df['Sample ID'].astype(str)
        # seq_report_df['isolate'] = "SARS-CoV-2/human/USA/" + seq_report_df["Sample ID"].astype(str) + "/" + pd.DatetimeIndex(seq_report_df['Test date']).strftime('%Y')
        return seq_report_df

    def _prep_barcode_df(self):
        """Creates df from `barcode_map`, if provided. Maps barcode to sample_name"""

        if self.barcode_map == None:
            return pd.DataFrame()
        bc = pd.read_csv(self.barcode_map,sep="\t",names=["barcode","sample_name"])
        if self.controls and type(self.controls)==str:
            bc = bc[bc["sample_name"].str.contains(self.controls)==False]
        return bc

    def _prep_gisaid_df(self):
        """Creates df from `gisaid_log`, if provided. Maps Virus name to accession"""

        # merge in gisaid accession details & update sample_name to include our location-based prefix
        if self.gisaid_log == None:
            return pd.DataFrame()
        else:
            gisaid_df = pd.read_csv(self.gisaid_log, sep=';', header=None, names=["gisaid_virus_name", "gisaid_accession"],usecols=range(2)).dropna()
            for col in gisaid_df.columns: # remove extra spaces in values
                gisaid_df[col] = gisaid_df[col].apply(lambda x: str(x).strip())
                gisaid_df["sample_name"] = gisaid_df["gisaid_virus_name"].apply(self.virus_to_sample_name)
                gisaid_df = gisaid_df[~gisaid_df["gisaid_virus_name"].str.contains("submissions")]
                gisaid_df = gisaid_df[~gisaid_df["sample_name"].isin(self._findExcludables())]
            return gisaid_df

    def _ensure_new_names_only(self,df):
        """Raises ValueError if any 'sample_name' in `df` was already-submitted"""

        if not self.submitted_samples_file:
            logging.warning(f"`submitted_samples` file not provided in `config_file`. Not checking for previously submitted samples.")
            return
        if not Path(self.submitted_samples_file).exists():
            logging.warning(f"Couldn't find file: {self.submitted_samples_file}. Not checking for previously submitted samples.")
            return
        submitted = pd.read_csv(self.submitted_samples_file,header=None,names=["already_submitted"])["already_submitted"].unique()
        attempted_duplicate = df[df["sample_name"].isin(submitted)]
        if len(attempted_duplicate) > 0:
            raise ValueError(f"The following sample(s) appear to have already been submitted.\n\n{attempted_duplicate}\n\n * To submit anyway,  remove that sample name from the file ({self.submitted_samples_file}) or temporarily comment out the variable `submitted_samples` in the `config_file`. \n * To exclude a sample from submission, add it to {self.exclude_file}")

    def _findExcludables(self):
        """Returns list of samples to exclude (if found in file `self.exclude_file`)"""

        excludables = []
        if self.exclude_file.exists():
            with self.exclude_file.open() as fh:
                for line in fh:
                    line = line.strip()
                    if line:
                        excludables.append(line)
        return excludables

    def _add_bioproject_if_needed(self,df:pd.DataFrame):
        """Adds default BioProject accession to field "bioproject_accession" if empty or absent.

        If a valid-looking accession exists, it will not be overwritten.
        If `test_dir==True`: sets all accessions to default test accession.

        Args:
            df (DatFrame): the dataframe to which to add accessions

        Raises:
            ValueError: if accessions appear invalid
        """

        default_production_accession = self.bioproject_presets.get("bioproject_accession")   # vvv a random test-dir bioproject accession
        default_test_accession = self.bioproject_presets.get("bioproject_test_accession","PRJNA553747")
        if self.test_dir:
            if not isBioProjectAccession(default_test_accession):
                raise ValueError(f"Default test BioProject accession '{default_test_accession}' in config_file '{self.config_file}' does not fit the standard format (PRJ[D|E|N]xxxxxx or PSUBxxxxxx).")
            df["bioproject_accession"] = default_test_accession
        else:
            if not isBioProjectAccession(default_production_accession):
                raise ValueError(f"Default BioProject accession '{default_production_accession}' in config_file '{self.config_file}' does not fit the standard format (PRJ[D|E|N]xxxxxx or PSUBxxxxxx).")
            if not "bioproject_accession" in df.columns:
                df["bioproject_accession"] = default_production_accession
            else:
                df = df.fillna(default_production_accession)
                # ensure all accessions look like accessions
                checkBioProjectAccessions(df["bioproject_accession"])
        return df

    def _prep_biosample_df(self,ignore_dates=False):
        """Creates BioSample df based on data in `seq_report` and `config_file`"""

        # act as a getter, if exists
        if type(self.biosample["df"])==pd.DataFrame:
            return self.biosample["df"]
        # create new df
        if self.use_existing:
            logging.info("Reading in sra df - not creating it fresh")
            biosample_df = pd.read_csv(self.biosample["tsv"],sep="\t")
            actual_cols = list(biosample_df.columns)
        else:
            # ^^ return existing df or vv create anew
            seq_report_df = self.seq_report_df
            # get variables from seq_report_df and set column headers to match those in template
            starting_cols = [x.strip() for x in self.biosample_presets.get("all_cols","").split(",") if x.strip() in seq_report_df.columns]
            if not starting_cols:
                starting_cols = [col for col in ["Sample #","Ct N gene","Test date","Sample ID","bioproject_accession","Primer Scheme"] if col in seq_report_df.columns]
                biosample_df = seq_report_df[starting_cols].rename(columns={"Test date":"collection_date","Sample #":"sample_name","Sample ID":"sample_title","Ct N gene":"sars_cov_2_diag_pcr_ct_value_1"})
            else:
                biosample_df = seq_report_df

            # add derived attributes
            # ensure these are added to the df and list of expected columns
            ct_cols = []
            for i,col in self.sars_cov_2_diag_pcr_ct_values.items():
                biosample_df[f'sars_cov_2_diag_pcr_ct_value_{i}'] = seq_report_df[col]
                ct_cols.extend([f'sars_cov_2_diag_pcr_ct_value_{i}',f'sars_cov_2_diag_gene_name_{i}'])

            # merge in gisaid accessions
            # biosample["sample_name_short"] = biosample['sample_name'].astype(str).apply(lambda x: x.split("-")[-1])
            # biosample = pd.merge(biosample,gisaid,on="sample_name_short",how='outer')
            if not self.gisaid_df.empty:
                # biosample = pd.merge(biosample,gisaid,on="sample_name_short",how='right')
                biosample_df = pd.merge(biosample_df,self.gisaid_df,on="sample_name",how='right')
                biosample_df = biosample_df[biosample_df["gisaid_accession"].notna()]
            else:
                biosample_df['gisaid_accession'] = 'missing'
            if ignore_dates:
                biosample_df['isolate'] = ''
            else:
                biosample_df['isolate'] = "SARS-CoV-2/human/USA/" + biosample_df["sample_name"].astype(str) + "/" + pd.DatetimeIndex(biosample_df['collection_date']).strftime('%Y')
            biosample_df["collection_date"] = pd.to_datetime(biosample_df["collection_date"]) # ensures dates are in the desired format
            biosample_df['gisaid_accession'] = biosample_df['gisaid_accession'].apply(lambda x: x if str(x)!='nan' else 'missing')

            # add config defaults/overrides
            for k,col in self.biosample_presets.items():
                if not k in ("bioproject_accession","bioproject_test_accession"):
                    biosample_df[k] = col

            # use test bioproject for all samples (if --test_dir==True)
            # use default bioproject for any samples that don't have one (if --test_dir==False)
            biosample_df = self._add_bioproject_if_needed(biosample_df)

            # finalize columns
            # all_bs_attr_ordered = ["sample_name","sample_title","bioproject_accession","organism","collected_by","collection_date","geo_loc_name","host","host_disease","isolate","isolation_source","antiviral_treatment_agent","collection_device","collection_method","date_of_prior_antiviral_treat","date_of_prior_sars_cov_2_infection","date_of_sars_cov_2_vaccination","exposure_event","geo_loc_exposure","gisaid_accession","gisaid_virus_name","host_age","host_anatomical_material","host_anatomical_part","host_body_product","host_disease_outcome","host_health_state","host_recent_travel_loc","host_recent_travel_return_date","host_sex","host_specimen_voucher","host_subject_id","lat_lon","passage_method","passage_number","prior_sars_cov_2_antiviral_treat","prior_sars_cov_2_infection","prior_sars_cov_2_vaccination","purpose_of_sampling","purpose_of_sequencing","sars_cov_2_diag_gene_name_1","sars_cov_2_diag_gene_name_2","sars_cov_2_diag_gene_name_3","sars_cov_2_diag_pcr_ct_value_1","sars_cov_2_diag_pcr_ct_value_2","sars_cov_2_diag_pcr_ct_value_3","sequenced_by","vaccine_received","virus_isolate_of_prior_infection","description"]
            # set attributes by package (for SARS-CoV-2 packages, only)
            if self.biosample_presets["package"] == "SARS-CoV-2.cl.1.0":
                required_bs_attr = ["sample_name","organism","collected_by","collection_date","geo_loc_name","host","host_disease","isolate","isolation_source"]
            elif self.biosample_presets["package"] == "SARS-CoV-2.wwsurv.1.0":
                required_bs_attr = ["sample_name","organism","collected_by","collection_date","geo_loc_name","isolation_source"]
            else: # This assumes the user hasn't forgotten any necessary columns. NCBI submission will fail if important attributes are missing, anyway...
                required_bs_attr = []
            all_bs_cols = [x.strip() for x in self.biosample_presets.get("all_cols","").split(",") if x]
            if not all_bs_cols:
                all_bs_cols = ['sample_name','sample_title','bioproject_accession','organism','collected_by','collection_date','geo_loc_name','host','host_disease','isolate','isolation_source','antiviral_treatment_agent','collection_device','collection_method','date_of_prior_antiviral_treat','date_of_prior_sars_cov_2_infection','date_of_sars_cov_2_vaccination','exposure_event','geo_loc_exposure','gisaid_accession','gisaid_virus_name','host_age','host_anatomical_material','host_anatomical_part','host_body_product','host_disease_outcome','host_health_state','host_recent_travel_loc','host_recent_travel_return_date','host_sex','host_specimen_voucher','host_subject_id','lat_lon','passage_method','passage_number','prior_sars_cov_2_antiviral_treat','prior_sars_cov_2_infection','prior_sars_cov_2_vaccination','purpose_of_sampling','purpose_of_sequencing']
                all_bs_cols.extend(ct_cols)
                all_bs_cols.extend(['sequenced_by','vaccine_received','virus_isolate_of_prior_infection','description'])
            if ignore_dates:
                required_bs_attr = remove_item("collection_date",required_bs_attr)
                all_bs_cols = remove_item("collection_date",all_bs_cols)
            actual_cols = [col for col in all_bs_cols if col in biosample_df.columns]
            check_missing(actual_cols,required_cols=required_bs_attr,name="BioSample")
            # biosample_df = biosample_df[actual_cols] # NOTE doing this later (before writing tsv and xml)
            ensureAllSamplesHaveNames(biosample_df)
            # filter out unwanted samples, if given
            biosample_df = biosample_df[~biosample_df['sample_name'].isin(self._findExcludables())]
            # ensure length of biosample matches length of barcodes available (excluding controls)
            if not self.barcode_df.empty:
                barcode_df_filtered = self.barcode_df[~self.barcode_df['sample_name'].isin(self._findExcludables())]
                extras = set(barcode_df_filtered["sample_name"].unique()) - set(biosample_df["sample_name"].unique())
                if len(extras) > 0:
                    logging.warning(dedent(f"""
                        Samples exist in barcode file that aren't found in BioSample file:
                        \t{extras}
                        Any sample can be skipped by writing it in a line by itself in a file called:
                        {self.exclude_file}
                        If {'these samples' if len(extras)!=1 else 'this sample'} should be excluded from submissions, you can use this command:"""))
                    for sample in extras:
                        print(f"echo '{sample}' >> '{self.exclude_file}'")
                    warn("")
            # ensure all samples present have not been previously submitted
            self._ensure_new_names_only(biosample_df)

        logging.info("BioSample")
        logging.info(biosample_df)
        self.biosample["final_cols"] = actual_cols

        self.biosample["df"] = biosample_df
        return biosample_df

    # def seekAccessionsDict(self):
    def get_accessions(self,as_dict=False,as_df=False,biosample_accessions:Path=None,require_biosample=True):
        """Returns dict or df mapping sequence id to biosample accession

        Args:
            as_dict (bool, optional): _description_. Defaults to False.
            as_df (bool, optional): _description_. Defaults to False.
            biosample_accessions (Path, optional): path to accession file. Defaults to None.
                If provided, uses file to determine accessions
                If not, recursively searches all report*.xml files in `report_dir` for one containing BioSample accessions

        Raises:
            ValueError: if requesting two different output types
            FileNotFoundError: if file provided but doesn't exist
        """

        if as_dict and as_df: raise ValueError("Choose only one of `as_dict` or `as_df`")

        # read accessions from spreadsheet
        if biosample_accessions!=None:
            df = pd.read_csv(biosample_accessions,sep="\t")
            df = df.rename(columns={
                "Accession":"BioSample","accession":"BioSample",
                "Sample name":"Sequence_ID","sample_name":"Sequence_ID"})
            if not "BioSample" in df.columns:
                warn(f"Unexpected headers in accession file:\n\t{biosample_accessions}\nThe fields 'accession' and 'sample_name' must be present.")
            if as_dict:
                return dict(zip(df["Sequence_ID"],df["BioSample"]))
            elif as_df:
                return df
            
        # search report*.xml files for one containing BioSample accessions
        else:
            accessions = {}
            # verify reports exist
            if not Path(self.report_dir).exists():
                if require_biosample:
                    raise FileNotFoundError(f"BioSample report directory not found:\n{self.report_dir}")
            else:
                # search for report with accessions
                for file in sorted([f for f in self.report_dir.glob("*/report.*") if f.name != "report.xml"],reverse=True):
                    logging.info(f"Checking {file} for accessions")
                    report = Report(file)
                    if report.biosamplesOk():
                        accessions_in_file = report.getAccessionDict(by_sample_name=True)
                        if accessions_in_file:
                            logging.info("Found accessions in",file)
                            accessions.update(accessions_in_file)
            if as_df:
                return pd.DataFrame(accessions.items(),columns=["Sequence_ID","BioSample"])
            elif as_dict:
                return accessions

    def add_biosample_2_genbank(self,biosample_accessions=None,require_biosample=False):
        """Adds BioSample accessions to GenBank TSV
        
        Args:
            biosample_accessions (str|Path): Path to tsv mapping biosample accessions to sample_name
            require_biosample (bool, optional): If True, raises AttributError if accessions can't be added.
        """

        # read in dfs
        if self.use_existing:
            logging.warning("Reading in genbank df - not creating it fresh")
            genbank_df = pd.read_csv(self.genbank["tsv"],sep="\t")
        else:
            genbank_df = self.genbank["df"]
            if type(genbank_df) != pd.DataFrame:
                self.prep_dfs(add_biosample=False)
        
        # only add biosamples if not already present
        if genbank_df["BioSample"].isnull().values.any() or "Missing" in genbank_df["BioSample"].unique():
            g_cols = genbank_df.columns
            genbank_df = genbank_df.drop(columns="BioSample")
            logging.info("reading in accessions")
            accessions_df = self.get_accessions(as_df=True,biosample_accessions=biosample_accessions,require_biosample=require_biosample)

            logging.info("merging in accessions")
            genbank_df = genbank_df.merge(accessions_df,on="Sequence_ID",how="outer")[g_cols]
            genbank_df = genbank_df[genbank_df['note'].notna()]

        if require_biosample:
            # ensure BioSamples were actually added to all samples
            if genbank_df["BioSample"].isnull().values.any() or "Missing" in genbank_df["BioSample"].unique():
                raise AttributeError(f"BioSample accessions could not be added. They may not yet exist in logs. Try running with flags `ftp --check` or `--biosample_accessions`")

        logging.info(genbank_df)
        self.genbank["df"] = genbank_df
        return genbank_df

    def get_gisaid_submitted_samples(self):
        """Returns list of samples that were submitted to gisaid (based on `gisaid_log`) or "all" """

        if self.gisaid_submitted_samples == None:
            if not self.gisaid_df.empty:
                return list(self.gisaid_df["sample_name"])
            else: return "all"
        else: return self.gisaid_submitted_samples

    def _prep_genbank_df(self,add_biosample=False,biosample_accessions=None,ignore_dates=False,require_biosample=False):
        """Creates GenBank df or adds BioSample accessions to df from existing GenBank TSV

        Args:
            add_biosample (bool, optional): A flag to locate and add BioSample accessions to existing GenBank TSV. Defaults to False.
            biosample_accessions (str | Path, optional): path to a file downloaded from NCBI containing BioSample accessions
                if not provided, accessions will be sought from report.xml files
            ignore_dates (bool, optional): A flag to drop date columns from df
            require_biosample (bool, optional): If True, raises AttributError if accessions can't be added.
        """

        if self.use_existing:
            logging.warning("Reading in genbank df - not creating it fresh")
            genbank_df = pd.read_csv(self.genbank["tsv"],sep="\t")
        else:
            cols = ['sample_name','organism','geo_loc_name', 'host', 'isolate', 'collection_date', 'isolation_source', 'bioproject_accession', 'gisaid_accession']
            if ignore_dates: cols = remove_item('collection_date',cols)
            genbank_df: pd.DataFrame = self.biosample["df"].copy()[cols]
            genbank_df['biosample_accession'] = 'Missing'
            genbank_df['gisaid_accession'] = genbank_df['gisaid_accession'].apply(lambda x: 'GISAID accession: ' + str(x).strip() if str(x)!='missing' else "")
            submittable = self.get_gisaid_submitted_samples()
            # filter to gisaid-only
            if not submittable=="all":
                genbank_df = genbank_df[genbank_df["sample_name"].isin(submittable)]
            new_cols = cols + ['biosample_accession']
            genbank_df = genbank_df[new_cols]
            genbank_df = genbank_df.rename(columns={'sample_name':'Sequence_ID','geo_loc_name':'country', 'host':'host', 'isolate':'isolate', 'collection_date':'collection-date', 'isolation_source':'isolation-source', 'biosample_accession':'BioSample', 'bioproject_accession':'BioProject', 'gisaid_accession':'note'})
        
        logging.info("GENBANK")
        logging.info(genbank_df)
        self.genbank["df"] = genbank_df
        self.genbank["final_cols"] = list(genbank_df.columns)

        if add_biosample:
            genbank_df = self.add_biosample_2_genbank(biosample_accessions,require_biosample=require_biosample)
            return
        
        return genbank_df
    
    def _offer_skip_option(self,samples):
        """Guides user on how to mark desired samples to be excluded from submission.

        Args:
            samples (Collection): Samples that are missing from one input spreadsheet or another
        """

        output = ["Any samples that should not be submitted can be written on a line by themselves in a file called",
            f"{self.exclude_file}",
            "You can use this command:"]
        for sample_name in samples:
            output += [f'echo "{sample_name}" >> {self.exclude_file}']
        return "\n".join(output)

    def _get_fastq_file(self,sample_name):
        """Verifies and returns single fastq file basename for SRA
        
        Args:
            sample_name (str): the sample for which to seek fastq

        Raises:
            FileNotFoundError: if single fastq file not found for `sample_name`
        """

        possible_files = list(self.fastq_dir.glob(f"*{sample_name}*.fastq"))
        if len(possible_files) == 1:
            file_path:Path = possible_files[0]
        elif len(possible_files) > 1:
            raise FileNotFoundError(f"Too many potential fastq files for sample '{sample_name}' in {self.fastq_dir}")
        else:
            raise FileNotFoundError(f"Can't find expected fastq file for {sample_name} in {self.fastq_dir}\n{self._offer_skip_option([sample_name])}")
        if file_path.exists():
            return file_path.name
        else:
            raise FileNotFoundError(f"expected fastq file {file_path} does not exist.\n{self._offer_skip_option([sample_name])}")
            # return None

    def _get_paired_fastq_files(self,sample_name):
        """Returns names of paired fastq files for SRA
        
        Args:
            sample_name (str): the sample for which to seek fastq

        Raises:
            FileNotFoundError: if paired fastq files not found for `sample_name`
        """

        possible_files:list[Path] = list(self.fastq_dir.glob(f"*{sample_name}*"))
        if len(possible_files) == 2:
            if possible_files[0].name.endswith("fastq") and possible_files[1].name.endswith("fastq"):
                return possible_files
        elif len(possible_files) == 0:
            raise FileNotFoundError(f"Can't locate fastqs in {self.fastq_dir} via sample name '{sample_name}'.\n{self._offer_skip_option(sample_name)}")
        elif len(possible_files) == 1:
            if possible_files[0].name.endswith("fastq"):
                raise FileNotFoundError(f"Expected directory containing paired fastq files but found single fastq:\n\t'{possible_files[0]}'")
            poss_directory = possible_files[0]
            possible_files = list(poss_directory.glob("*.fastq"))
            if len(possible_files) == 2:
                if possible_files[0].name.endswith("fastq") and possible_files[1].name.endswith("fastq"):
                    return possible_files
            else:
                raise FileNotFoundError(f"Can't locate any paired '*.fastq' files in {poss_directory}")

    def _add_primer_schemes(self,df):
        """Adds primer scheme column to df if possible/needed
        
        Args:
            df (DataFrame): the dataset to which to add Primer Scheme info

        Raises:
            AttributeError: if missing required fields/attributes
            ValueError: if provided scheme is not in `allowed_schemes`
        """

        if "Primer Scheme" in df.columns:
            pass
        # attempt to get scheme from primer map file
        elif self.primer_map:
            primers = pd.read_csv(self.primer_map).rename(
                columns={
                "Sample #":"sample_name","Seq ID":"sample_name",
                "primer_scheme":"Primer Scheme",
                "primer scheme":"Primer Scheme",
                })
            if "Primer Scheme" in primers.columns:
                primers = primers[["sample_name","Primer Scheme"]]
                df = df.merge(primers,on="sample_name")
                try_others = False
                primer_map_failed = False
            else: # if Primer Scheme column missing from file
                if not self.primer_scheme:
                    raise AttributeError(f"'Primer Scheme' must be a field in the `primer_map` file '{self.primer_map}'. Alternatively, provide `primer_scheme` as an argument.")
                else:
                    try_others = True
                    primer_map_failed = True
        if try_others == True:
            # if no primer_map or 'Primer Scheme' column isn't found therein, use default for all samples
            if self.primer_scheme:
                df["Primer Scheme"] = self.primer_scheme
            # If there's only one possible scheme in config, use it
            elif len(self.allowed_schemes) == 1:
                df["Primer Scheme"] = self.allowed_schemes[0]
            else: raise AttributeError(f"'Primer Scheme' must be a field in your `seq_report` file '{self.primer_scheme}' \n or else one of the arguments `--primer_map` or `--primer_scheme` must be provided")
        
        # verify primer schemes are allowed based on config
        if len(self.allowed_schemes) == 0: raise AttributeError(f"'allowed_schemes' must be specified in `config_file` '{self.config_file}'")
        extra_schemes = set(df["Primer Scheme"].unique()) - set(self.allowed_schemes)
        if len(extra_schemes) > 0:
            raise ValueError(f"The following scheme(s) are not among the `allowed_schemes` listed in the `config_file` '{self.config_file}':\n{extra_schemes}")

        return df

    def _addFilenames(self,df:pd.DataFrame):
        """Adds paths to filenames for each sample (pair reads get extra filename column)
        
        Args:
            df (DataFrame): dataframe to which to add filenames
        """

        if "illumina" in self.sra_presets["platform"].lower():
            df[["filename","filename2"]] = df['sample_name'].apply(self._get_paired_fastq_files,result_type="expand")
        else:
            df["filename"] = df["sample_name"].apply(self._get_fastq_file)
        return df

    def _prep_sra_df(self):
        """Returns DataFrame of SRA attributes derived from BioSample data or retrieved from `config_file`"""

        # act as a getter, if exists
        if type(self.sra["df"])==pd.DataFrame:
            return self.sra["df"]
        # create new df
        if self.use_existing:
            logging.warning("Reading in sra df - not creating it fresh")
            sra_df = pd.read_csv(self.sra["tsv"],sep="\t")
        else:
            # get repeated fields from BioSample
            sra_df = self._prep_biosample_df().copy()
            cols = ['bioproject_accession','sample_name','sample_title']
            if "Primer Scheme" in sra_df.columns:
                cols.append("Primer Scheme")
            sra_df = sra_df[cols]
            sra_df = sra_df[~sra_df['sample_name'].isin(self.controls.split("|"))] # weed out controls - we don't want to submit those
            submittable = self.get_gisaid_submitted_samples()
            if not self.barcode_df.empty:
                sra_df = sra_df.merge(self.barcode_df,on='sample_name',how="outer")
            if not submittable=="all":
                sra_df = sra_df[sra_df["sample_name"].isin(submittable)]
            else:
                sra_df = sra_df[sra_df["sample_name"].notna()]

            # add derived attributes
            sra_df["library_ID"] = sra_df["sample_name"]
            # remove any sampes that are supposed to be excluded
            if self.exclude_file.exists():
                with self.exclude_file.open() as fh:
                    sra_df = sra_df[~sra_df["sample_name"].isin(set([line.strip() for line in fh]))]
            sra_df = self._addFilenames(sra_df)
            sra_df = sra_df.dropna(subset=["filename"])
            sra_df = sra_df[sra_df["filename"]!=None]
            sra_df = self._add_primer_schemes(sra_df)
            for col in ("amplicon_PCR_primer_scheme","design_description","sequencing_protocol_name"):
                sra_df[col] = sra_df["Primer Scheme"].apply(lambda scheme: self.protocol_scheme[scheme][col])
            sra_df = sra_df.drop(columns=["Primer Scheme"])

            # add `config_file` defaults/overrides
            for k,v in self.sra_presets.items():
                sra_df[k] = v

        # finalize columns
        required_sra_cols = ["bioproject_accession","sample_name","library_ID","title","library_strategy","library_source","library_selection","library_layout","platform","instrument_model","design_description","filetype","filename"]
        all_sra_cols = ['bioproject_accession','sample_name','library_ID','title','library_strategy','library_source','library_selection','library_layout','platform','instrument_model','design_description','filetype','filename','filename2','amplicon_PCR_primer_scheme','amplicon_size','sequencing_protocol_name','raw_sequence_data_processing_method','dehosting_method','sequence_submitter_contact_email']
        extra_cols = [col for col in sra_df.keys() if col not in set(all_sra_cols)]
        all_sra_cols = all_sra_cols
        # actual_cols = [col for col in all_sra_cols+['sample_name_short'] if col in sra.columns] + extra_cols
        actual_cols = [col for col in all_sra_cols if col in sra_df.columns] + extra_cols

        # NOTE: this only allows for up to two filenames, at the moment
        # sra_df = sra_df[actual_cols] # NOTE doing this later (before writing tsv and xml)
        check_missing(actual_cols,required_sra_cols,"SRA")

        logging.info("SRA")
        logging.info(sra_df)
        self.sra["df"] = sra_df
        self.sra["final_cols"] = actual_cols
        return sra_df

    # def merge_dfs(self,data,df_names):
    def _prep_merged_df(self):
        """Returns combined DataFrame of all samples/attributes"""

        biosample_df, sra_df, genbank_df = [getattr(self,dataset)["df"] for dataset in ("biosample","sra","genbank")]
        biosample_df:pd.DataFrame
        sra_df:pd.DataFrame
        genbank_df:pd.DataFrame
        # create merged df (before genbank['sample_title'] gets dropped)
        s1=set(biosample_df.columns)
        s2=set(sra_df.columns)
        sra_extras = remove_item('sample_name',s1.intersection(s2))
        merged_df = biosample_df.merge(sra_df.drop(columns=sra_extras),on='sample_name',how='outer')
        if not genbank_df.empty:
            s3=set(genbank_df.columns)
            genbank_extras = remove_item('sample_name',s1.intersection(s3))
            merged_df = merged_df.merge(genbank_df.drop(columns=genbank_extras),left_on='sample_name',right_on='Sequence_ID',how='outer')
        self.merged["df"] = merged_df
        return merged_df

    def check_existing_tsvs(self):
        """Raises FileNotFoundError if TSVs don't already exist"""

        if self.use_existing:
            for dataset in (self.sra, self.biosample, self.genbank):
                if not dataset["tsv"].exists():
                    raise FileNotFoundError(f"At least one `file_prep` TSV not found. Try again without the `use_existing` flag.")

    def prep_dfs(self,add_biosample=False,ignore_dates=False,biosample_accessions=None,require_biosample=False):
        """Prepares spreadsheets for NCBI submission portal (https://submit.ncbi.nlm.nih.gov/)

        Args:
            sra (bool, optional): A flag to create SRA TSV. Defaults to True.
            biosample (bool, optional): A flag to create BioSample TSV. Defaults to True.
            genbank (bool, optional): A flag to create GenBank TSV. Defaults to True.
            add_biosample (bool, optional): A flag to add biosample accessions to GenBank TSV. Defaults to False.
            ignore_dates (bool, optional): A flag to drop date columns. Defaults to False.
            require_biosample (bool, optional): If True, raises AttributError if accessions can't be added.
        """

        self.check_existing_tsvs()
        self._prep_biosample_df(ignore_dates)
        self._prep_sra_df()
        self._prep_genbank_df(add_biosample=add_biosample,biosample_accessions=biosample_accessions,ignore_dates=ignore_dates,require_biosample=require_biosample)
        self._prep_merged_df()

    ###############  ^^^^  tsv_prep  ^^^^  ###############
    ###############  vvvv  xml_prep  vvvv  ###############

    ###############  ^^^^  xml_prep  ^^^^  ###############
    ###############  vvvv  zip_prep  vvvv  ###############

    def verifyUnsubmittable(self,samples_to_maybe_exclude,warning=""):
        """Prints warning and recommends marking samples for exclusion if needed"""

        if len(samples_to_maybe_exclude) != 0:
            logging.warning(warning)
            self._offer_skip_option(samples_to_maybe_exclude)
            exit(1)

    def getAllowedSamples(self,all_samples=None):
        """Returns collection of allowed samples from gisaid logfile or else all_samples
        
        Args:
            all_samples (Collection): default samples if no gisaid-submitted-samples exist to limit by
        """

        allowed_samples = self.get_gisaid_submitted_samples()
        if allowed_samples == "all":
            logging.warning(f"No `gisaid_log` provided, so assuming all samples from `seq_report` should be submitted excluding any found in {self.exclude_file}")
            # NOTE: filtering out samples from `self.exclude_file` has already been done in TSV prep
            return all_samples
        gisaid_extras = set(allowed_samples) - set(all_samples)
        gisaid_missing = set(all_samples) - set(allowed_samples)
        self.verifyUnsubmittable(gisaid_extras,"Sample(s) found in `gisaid_log` that lack metadata.")
        self.verifyUnsubmittable(gisaid_missing,"Metadata exists for sample(s) that can't be found in `gisaid_log`.")
        return set(allowed_samples)
    
    def generate_updated_fasta_records(self):
        """Yields fasta records with only first part of header (drops "/medaka*")"""

        p = Path(self.fasta)
        if not p.exists():
            raise FileNotFoundError(f"`fasta` file not found: {self.fasta}")
        for record in SeqIO.parse(self.fasta,"fasta"):
            # update id # TODO: ensure this splitting is universal
            record.id = record.id.split("/")[0]
            yield record

    def write_genbank_fasta(self,outHandle):
        """Writes GenBank fasta seqs to provided outHandle"""

        allowed_samples = self.getAllowedSamples(self.biosample["df"]["sample_name"].unique())
        added_samples = set()
        # write fasta
        for record in self.generate_updated_fasta_records():
            if record.id in allowed_samples:
                added_samples.add(record.id)
                outHandle.write(f">{record.id}\n")
                # trim leading/trailing Ns
                record.seq = record.seq.strip("N")
                outHandle.write(f"{record.seq}\n")
        # ensure all expecte samples were found and added to genbank fasta
        if len(allowed_samples) > len(added_samples):
            self.verifyUnsubmittable(allowed_samples-added_samples,f"Not all samples found in `fasta` file: {self.fasta}.")
    
    def write_genbank_fasta_to_zip(self,zfh):
        """Writes genbank fasta to genbank.zip/seqs.fsa"""
        
        with StringIO() as memFh:
            self.write_genbank_fasta(memFh)
            zfh.writestr("seqs.fsa", memFh.getvalue())

    def write_genbank_seqs_src(self,zfh,biosample_accessions=None):
        """Writes genbank.tsv to genbank.zip/seqs.scr

        Args:
            biosample_accessions (str|Path): Path to tsv mapping biosample accessions to sample_name
        """
        
        # add accessions to GenBank TSV
        self.add_biosample_2_genbank(biosample_accessions,require_biosample=True)
        genbank_tsv = self.outdir.joinpath("genbank-metadata.tsv")
        df_2_tsv(self.genbank["df"],genbank_tsv)
        # copy updated TSV to Zip file
        with genbank_tsv.open() as tsv:
            tsv_string = tsv.read()
            zfh.writestr("seqs.src", tsv_string)

    def write_genbank_template(self,zfh):
        """Writes genbank.tsv to genbank.zip/seqs.scr"""

        template = getattr(self,"template",None) # seek template path from `config_file`
        if not template:
            # TODO: generate template from `config_file`
            raise Exception("Path to `template` file must be provided in `config_file`.")
        tentaive_template = Path(template).resolve()
        if tentaive_template.exists():
            template = tentaive_template
        else:
            if self.config_file.parent.joinpath(template).exists():
                template = self.config_file.parent.joinpath(template)
            else: raise FileNotFoundError("Can't find template: {template}")
        with template.open() as sbt, StringIO() as memFh:
            for line in sbt:
                if "Submission Title:None" in line:
                    line = line.replace(":None",f":{self.plate} GenBank")
                memFh.write(line)
            zfh.writestr("seqs.sbt",memFh.getvalue())
            # zfh.writestr("seqs.sbt",sbt.read())
        # zfh.writestr("seqs.sbt", makeTemplate())

    def write_genbank_comment_file(self,zfh):
        """Writes structured comment to genbank.zip/comment.cmt"""

        if self.genbank_presets["comment"]["include_comment"]:
            # create comment df based on config presets
            # by now, accessions have been added to genbank_df
            comment_df:pd.DataFrame = self.genbank["df"][["BioSample","Sequence_ID"]].rename(columns={"BioSample":"BioSample Accession","Sequence_ID":"sample_name"})
            for field,default in self.genbank_presets["comment"].get("data",{}).items():
                comment_df[field] = default
            # write it out
            with StringIO() as memFh:
                comment_df.to_csv(memFh,sep="\t",index=False)
                zfh.writestr("comment.cmt", memFh.getvalue())

    def write_genbank_zip(self,biosample_accessions=None):
        """Writes zipfile for GenBank submission

        Args:
            biosample_accessions (str|Path): Path to tsv mapping biosample accessions to sample_name
        """

        with ZipFile(self.genbank['zip_file'], "w") as zfh:
            print("\tWriting fasta --> seqs.fsa")
            self.write_genbank_fasta_to_zip(zfh)
            print("\tWriting genbank-metadata.tsv --> seqs.src")
            self.write_genbank_seqs_src(zfh,biosample_accessions)
            print("\tWriting template --> template.sbt")
            self.write_genbank_template(zfh)
            print("\tWriting structured comment --> comment.cmt")
            self.write_genbank_comment_file(zfh)

    ###############  ^^^^  zip_prep  ^^^^  ###############

    def write_genbank_xml(self):
        """Writes out XML for GenBank submission"""

        xml_maker = GenBank_Submission(self)
        xml_maker.write_xml(self.genbank["xml_file"])

    def write_sra_biosample_xml(self):
        """Writes out XML for SRA and/or BioSample submission"""

        xml_maker = SRA_BioSample_Submission(self)
        xml_maker.write_xml(self.sra["xml_file"])

    def write_presubmission_metadata(self,sra_only=False,require_biosample=False):
        """Prepares TSV and XML files

        This should be the go-to method for preparing and writing out files for initial submission to BioSample and/or SRA

        Args:
            action (str: "file_prep" | "add_biosample"): 
                "file_prep": Creates initial TSV and XML files
                "add_biosample": adds biosample accessions from previous submission to genbank tsv
            keep_tsvs (bool, optional): If False, TSVs will be deleted. Defaults to True.
            keep_xmls (bool, optional): If False, XMLs will be deleted. Defaults to True.
            require_biosample (bool, optional): If True, raises AttributError if accessions can't be added.
        """

        ignore_dates = True if sra_only else False

        # prepare/write TSVs
        self.prep_dfs(ignore_dates=ignore_dates,require_biosample=require_biosample,add_biosample=True)

        # final_dfs = {} # NOTE: not necessary?
        for dataset in (self.sra, self.biosample, self.genbank):
            df_2_tsv(
                df=finalize_df(df=dataset["df"],final_cols=dataset["final_cols"]),
                outfile=dataset["tsv"],name=dataset["name"]
            )

        # prepare/write XML
        print(f'Writing sra/biosample xml:\n   {self.sra["xml_file"]}')
        self.write_sra_biosample_xml()

    def write_genbank_submission_zip(self,biosample_accessions=None):
        """Prepares Zip file for GenBank submission.
        
        Args:
            biosample_accessions (str|Path): Path to tsv mapping biosample accessions to sample_name

        Contents:
          * genbank.xml
          * genbank.zip (seqs.sbt, seqs.fsa, seqs.src, [comment.cmt])
        """

        print(f"Writing genbank xml:\n   {self.genbank['xml_file']}")
        self.write_genbank_xml()
        print(f"Writing genbank zip:\n   {self.genbank['zip_file']}")
        self.write_genbank_zip(biosample_accessions)


    ### --- ^^^ --- file prep --- ^^^ --- ###

    ### --- vvv --- ncbi interaction --- vvv --- ###

    # NOTE: TODO: add in replaceTestAccessions() before submitting
        # only needed if people copy files from test_dir to non-test_dir...

    def _set_default_subdir(self,subdir):
        """Determines name of subdir
        
        Args:
            subdir (str): desired subdirectory name (just the basename, not path)

        Raises:
            AttributeError: if subdir or plate not provided
        """

        if subdir: self.subdir = subdir
        elif self.plate: self.subdir = self.plate
        else: raise AttributeError("Either `subdir` or `plate` must be provided.")

        return self.subdir


    def _set_ftp_vars(self):
        """Determines the paths to useful locations for FTP submission/checking"""

        self.submit_area = "submit/Test" if self.test_dir else 'submit/Production'
        self.submit_dir = Path(self.submit_area) / self.subdir
        self.initial_submission_dirs = []
        self.ftp_file_sizes = {}

        return self.submit_area,self.submit_dir,self.initial_submission_dirs,self.ftp_file_sizes

    def ncbiConnect(self) -> FTP:
        """Returns an FTP object connected to NCBI in the `submit_area`"""

        for arg_name,var in {"username (ncbiUser)":self.ncbiUser,"password (ncbiPass)":self.ncbiPass}.items():
            if not var:
                raise Exception(f"{arg_name} must be supplied on the command line, in the `config_file`, or as an environmental variable.")
        self.ftp = FTP(self.host,self.ncbiUser,self.ncbiPass)

        # navigate to submit/Test or submit/Production
        self.ftp.cwd(self.submit_area)
        logging.info(f"CWD = submit_area: {self.ftp.pwd()}\n")

    def ensureValidDatabase(self,db):
        """Raises ValueError if `db` is invalid"""

        if db not in self.valid_dbs + [None]:
            raise ValueError(f"The `db` {db} is invalid. Acceptable `db` values: {self.valid_dbs}")

    def ncbi_interact(self,action,db=None,attempt_num=1,simple=False):
        """Submits files depending on desired `db`
        
        Args:
            db (str): NCBI database to which to submit. Options: ["sra","gb","bs_sra","bs",None]
        """

        self.ensureValidDatabase(db)

        # connect to NCBI via ftp and go to Test or Production dir
        self.ncbiConnect()

        # set list of any existing submissions
        self.initial_submission_dirs = sorted(self.ftp.nlst())

        # submit or check
        if action == "submit":
            self._do_submit(db,attempt_num)
        elif action == "check":
            self._do_check(db,attempt_num,simple)

        self.ftp.quit()


    def submit(self,db,attempt_num=1):
        """Submits files depending on desired `db`
        
        Args:
            db (str): NCBI database to which to submit. Options: ["sra","gb","bs_sra","bs"]
        """

        self.ncbi_interact("submit",db,attempt_num)


    def check(self,db=None,attempt_num=1,simple=False):
        """Checks on submission status. If `db` specified, only checks the one. Otherwise checks on any found.
        
        Args:
            db (str): NCBI database to which to submit. Options: ["sra","gb","bs_sra","bs"]
        """

        self.ncbi_interact("check",db,attempt_num,simple)

    def enterLocalDir(self,directory):
        """Makes `directory`, if needed, and enters it"""

        directory = Path(directory)
        directory.mkdir(parents="True",exist_ok="True")
        logging.info(f"Moving to '{directory}'")
        os.chdir(directory)
    
    def upload(self,file:Path,outfile=None):
        """Uploads `file` to submission dir as outfile"""

        if self.test_mode:
            print(f"Test mode - skipping upload of {file}")
        else:
            print(f"Uploading {file}")
            if not outfile: outfile = file
            with io.open(file,"rb") as fh:
                self.ftp.storbinary(f"STOR {outfile}", fh)
            self.uploaded_something = True

    def upload_if_not_there(self,file:Path):
        "Only uploads files with name/filesize combinations that don't already exist in submission dir"
        if Path != type(file) == str : file = Path(file)
        else: raise Exception(f"Unexpected file type for {file}: {type(file)}")
        if not file.exists():
            raise FileNotFoundError(f"Can't upload non-existent file: {file}")
        # if expected file exists with expected size, don't upload
        outfile = "submission.xml" if file.name.endswith(".xml") else file.name
        # allow small (non-fastqs (xmls)) to be updated whether size changed or not (but don't re-upload fastqs)
        if "fastq" in outfile and outfile in self.ftp_file_sizes.keys() and int(self.ftp_file_sizes[outfile]) == int(file.stat().st_size):
            print(f"Skipping upload of existing file: {file.name}")
            return
        self.upload(file,outfile)

    def mark_submitted(self,sample_name):
        """Appends `sample_name` to `submitted_samples` file if `submitted_samples` provided (and not `test_mode`)"""

        if self.test_mode:
            print(f"Would be adding {sample_name} to {self.submitted_samples_file}")
            return
        if self.submitted_samples_file:
            with open(self.submitted_samples_file,"a") as out:
                out.write(f"{sample_name}\n")

    def set_ftp_file_sizes(self):
        """Sets self.ftp_file_sizes dict of sizes of files in current ftp directory"""

        # list all files in current dir
        mlsd_size = [x for x in self.ftp.mlsd(facts=["size"])]
        self.ftp_file_sizes = {}
        logging.info("File sizes at NCBI:")
        for name,details in mlsd_size:
            # write details to stdout
            logging.info(name,details,"size" in details.keys())
            # Only keep files that have size info available
            if "size" in details.keys():
                self.ftp_file_sizes[name]=details["size"]

    def mark_submit_ready(self):
        """Creates submit.ready file, if anything was submitted"""

        # this file indicates to NCBI that the submission is ready to go
        if not self.test_mode and self.uploaded_something:
            print("Creating empty submit.ready file\n")
            with io.StringIO() as emptyFh:
                emptyFh.write(u"")
                self.ftp.storbinary('STOR submit.ready', emptyFh)

    def _do_submit(self,db,attempt_num=1):
        """Submits files depending on desired `db`
        
        Args:
            db (str): NCBI database to which to submit. Options: ["sra","gb","bs_sra","bs"]
        """

        # make submission directory if needed (and go there)
        print(f"existing directories in {self.submit_area}:")
        print(self.initial_submission_dirs)
        subdir = self.setSubdir(db,attempt_num)
        if not subdir in self.initial_submission_dirs:
            print(f"\nMaking subdir: {subdir}")
            self.ftp.mkd(subdir)
        else:
            logging.warning(f"Subdirectory '{subdir}' already exists. "
                "We assume you are finishing an incomplete submission or updating a previous submission.")
            
        # move into submission directory
        self.ftp.cwd(subdir)
        logging.info(f"\nCWD = subdir: {self.ftp.pwd()}\n")

        # move to directory of files to upload
        logging.info("Moving to local dir:",self.outdir)
        self.enterLocalDir(self.outdir)
        logging.info("Current local files:")
        logging.info(os.listdir())

        # get file sizes as dict {filename:size}
        self.set_ftp_file_sizes()

        print("\nUploading files...\n")

        # upload specific files from ncbi directory (`--outdir`)
        # genbank
        if db == "gb":
            for file in ("genbank.zip","genbank.xml"):
                self.upload_if_not_there(file)
        
        # biosample/sra
        if db == "bs_sra":
            self.upload_if_not_there("sra_biosample.xml")
            # move to directory containing fastqs files to upload
            os.chdir(self.fastq_dir)
            # find and upload all (fastq) files in any column labeled "filename*"
            sra_df = self._prep_sra_df()

            fn_cols = [col for col in sra_df.columns if col.startswith("filename")]
            for i,row in sra_df.iterrows():
                sample_name = row["sample_name"]
                for col in fn_cols:
                    file = row[col]
                    self.upload_if_not_there(file)
                self.mark_submitted(sample_name)

        self.mark_submit_ready()

        # show current files at NCBI
        print(f"Current files in remote dir '{self.ftp.pwd()}':")
        print(self.ftp.nlst())

    ### ^^ submit
    ### vv check

    def getFileFact(self,fact,contains_str=None):
        """Returns dict of {filename:fact}
        
        Args:
            fact (str): detail to get. Options: [modify, perm, size, type, unique, unix.group, unix.mode, unix.owner]
            contains_str (str): if provided, limits to files with `contains_str` in their name
        """

        # get all filenames and desired fact
        # mlsd_tuples = [x for x in ftp.mlsd(facts=[fact])]
        # convert to dict {filename:mod}
        fact_dict = {}
        for name,details in self.ftp.mlsd(facts=[fact]):
            if type(contains_str)!=type(None) and not contains_str in name:
                continue # if limiting the list, only keep if name contains this string
            if not fact in details.keys():
                continue # only keep if fact exists for name
            # save file:fact pair
            if fact == "modify":
                fact_dict[name] = asDate(details[fact])
            else:
                fact_dict[name] = details[fact]
        return fact_dict

    def getReportModTimes(self,as_str=True):
        """Yields time-stamps for each report file in NCBI
        
        Args:
            as_str (bool):
                `True`: yields str: "filename   time_stamp"
                `False`: yields tuple (filename, time_stamp)
        """

        # this ignores the simlink file report.xml
        mod_times = self.getFileFact("modify","report.")
        logging.info("mod_times:",mod_times)
        for report_file in sorted(list(mod_times.keys())):
            if as_str: yield f"{report_file}\t{mod_times[report_file]}"
            else: yield (report_file,mod_times[report_file])

    def isRequestedAttempt(self,submission_dir,subdir,attempt_num=1): #TODO: deprecate?
        """Returns True if `submission_dir` is the one requested based on the args `subdir` and `attempt_num`
        
        No longer used.
        Args:
            submission_dir (str): must just be the basename of the file (not a full path)
            subdir (str): spedific subdirectory name of this submission
            attempt_num (int|str): the attempt number of this submission
        """

        expected_end = f"_{attempt_num}" if attempt_num > 1 else ""
        expected_ncbi_name = f"{subdir}{expected_end}"
        return submission_dir == expected_ncbi_name

    def enterNcbiSubdir(self,ncbi_dir):
        """Moves to remote NCBI submission directory
        
        Args:
            ncbi_dir (str): name of a subdirectory of interest within the current `submit_area`
        """

        logging.info(f"Moving from {self.ftp.pwd()} --> {self.submit_area}/{ncbi_dir}")
        self.ftp.cwd(ncbi_dir)
        return True

    def numReports(self):
        """Returns number of NCBI reports available in current submission dir"""

        return len(set(f for f in self.ftp.nlst() if f.startswith("report.")))

    def relativePathToPrevDir(self,prev_dir,current_dir):
        """Returns relative path to previous directory. Full paths must be used."""

        prev_dir,current_dir = str(prev_dir),str(current_dir)
        logging.info("Moving from",prev_dir,"-->",current_dir)
        # if going from a/b/c/d --> a/b
        if current_dir.startswith(prev_dir):
            remainder = current_dir[len(prev_dir):]
            num_pieces = len([x for x in remainder if x=="/"]) + 1
            return "/".join(num_pieces*[".."])
        elif prev_dir.startswith(current_dir):
            remainder = prev_dir[len(current_dir):].strip("/")
            return f"./{remainder}"
        else: raise Exception("More work needed here to determine relative paths")

    def downloadReportFiles(self,report_dir):
        """Downloads all report*.xml files"""

        # enter localdownload directory
        starting_dir = os.getcwd()
        self.enterLocalDir(report_dir)
        for file in self.ftp.nlst():
            if file.startswith("report.") and file.endswith(".xml"):
                print(f"Downloading '{file}' to '{os.getcwd()}/{file}'")
                self.ftp.retrbinary("RETR " + file, open(file, 'wb').write)
        # return to original dir
        logging.info(f"Moving to '{starting_dir}'")
        os.chdir(starting_dir)

    def isReportFile(self,file) -> bool:
        """Returns True for files that are named like report files (report*.xml)"""

        return file.startswith("report.") and file.endswith(".xml")

    def submissionFiles(self):
        """Returns list of all files excluding report files"""

        return [f for f in self.ftp.nlst() if not self.isReportFile(f)]

    def fullDatabaseReport(self):
        """Prints info about current files in NCBI submission dir and report modification times"""

        print("\n\t\tModification times:")
        for report_time in self.getReportModTimes():
            print(f"\t\t\t{report_time}")
        print("\n\t\tFiles (excluding report*xml files)")
        for f in self.submissionFiles():
            print(f"\t\t\t{f}")

    def fullStatusReport(self,report_file):
        """Prints comprehensive details from `report_file`"""

        report = Report(report_file)
        print(f"\n\tReport file: '{report_file}'")
        print(f"\n\t\tOverall status:\n\t\t\t{report.simpleReport()}")
        report_list = report.statusReport(self.test_dir)
        if len(report_list) == 0:
            print("\t\tNo more details yet available.")
        print("\n\t\tStatus breakdown:")
        for line in report_list:
            if line.endswith(" status report"): print() # add empty line before every db-specific report
            print(f"\t\t\t{line}")

    def setSubdir(self,db,attempt_num=1):
        """Ensures submission directory is correct for the given database if checking both"""

        logging.info(f"setting subdir based on subdir='{self.subdir}' db='{db}'")
        db_str = f"_{db}" if db else ""
        num = f"_{attempt_num}" if attempt_num > 1 else ""
        actual_subdir = f"{self.subdir}{db_str}{num}"
        return actual_subdir

    def analyze_and_report(self,db,relevant_submission_dirs,local_db_report_dir,attempt_num,subdir):
        """Seeks useful details about submission(s) if extant and prints findings"""

        if len(relevant_submission_dirs) == 0:
            print("No submissions found for the database",db)
            return
        elif len(relevant_submission_dirs) > 0:
            print(f"\nChecking on '{db}' submission")
        elif len(relevant_submission_dirs) > 1:
            print(f"\nChecking all related submission dirs:\n\t{relevant_submission_dirs}")

        # seek info about each submission related to this same submission (subdir)
        for submission_dir in relevant_submission_dirs:

            starting_dir = os.getcwd()
            self.enterNcbiSubdir(submission_dir)
            isEmpty = len(self.ftp.nlst()) == 0

            # if directory is empty, skip this stuff
            if not isEmpty: # (at miniumum, some submission files exist; maybe no reports, though)
                report_dir = local_db_report_dir / submission_dir # local dir for report downloads

                # download reports if any exist (warning - this will overwrite any in local_db_report_dir with same name as those being downloaded)
                if self.numReports() > 0:
                    self.downloadReportFiles(report_dir)
                else:
                    print(f"No submission reports available for '{submission_dir}'")
                    return

                # print overview of db submission attempt
                print(f"\n\t{submission_dir} - overview of submission attempt")
                self.fullDatabaseReport()

                # give full report for each report.*.xml file
                for report_file in sorted([f for f in report_dir.glob("report*.xml") if f.name!="report.xml"]):
                    self.fullStatusReport(report_file)

            elif isEmpty:
                print("No data yet submitted for",submission_dir)
            # back out of directory so we can move into the next one with relative path
            self.ftp.cwd("..")
            # return to original local dir
            os.chdir(starting_dir)

    def analyze_if_possible(self,db,attempt_num=1,simple=False):
        """Seeks useful details about submission and prints findings"""

        # ensure correct `subdir` name
        db_subdir = self.setSubdir(db,attempt_num)
        local_db_report_dir:Path = self.outdir / f"{db}_reports"
        # will only check submission directories related to this `subdir`
        relevant_submission_dirs = [dir for dir in self.initial_submission_dirs if dir.startswith(db_subdir)]
        if relevant_submission_dirs:
            logging.info("subdir:",db_subdir)
            logging.info("relevant_submission_dirs:",relevant_submission_dirs)
        if simple: # wean down to the one current submission directory (based on `attempt_num`)
            expected_end = str(attempt_num) if attempt_num > 1 else db_subdir
            relevant_submission_dirs = [dir for dir in relevant_submission_dirs if dir.endswith(expected_end)]

        self.analyze_and_report(db,relevant_submission_dirs,local_db_report_dir,attempt_num,db_subdir)

    def _do_check(self,db,attempt_num=1,simple=False):
        """Checks on submission status. If `db` specified, only checks the one. Otherwise checks on any found.
        
        Args:
            db (str): NCBI database to which to submit. Options: ["sra","gb","bs_sra","bs"]
        """

        if db == None:
            if simple:
                raise AttributeError(f"If `simple` is specified, `db` must also be provided")
            else:
                # TODO: remove None as an option
                dbs = self.valid_dbs + [None]
        else: dbs = [db]

        # check submissions for anyspecified db
        for db in dbs:
            self.analyze_if_possible(db,attempt_num,simple)

    ### --- ^^^ --- ncbi interaction --- ^^^ --- ###
