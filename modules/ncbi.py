#!/usr/bin/env python3
"""A main class for preparing data and interacting with ncbi.
"""

import datetime
from io import StringIO
import os
import textwrap
from zipfile import ZipFile
import pandas as pd
from modules.helpers import *
from modules.report import Report
from modules.xml_format import SRA_BioSample_Submission,GenBank_Submission

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
# self.simple = simple
# self.attempt_num = attempt_num

class NCBI:
    def __init__(self,plate,fastq_dir,seq_report,outdir,barcode_map=None,
                 config=None,controls=None,gisaid_log=None,primer_map=None,
                 primer_scheme=None,fasta=None,ncbiUser=None,ncbiPass=None,
                 host=None,test_dir=False,test_mode=False) -> None:
        """Creates an object that can prepare, submit, and track csv/xml submission files 

        Args:
            config: (str, Path, None) Location of the config file to use
        """
        # get any variables from config
        cf,self.config_file = getConfig(config)
        self.template = cf["template"]
        self.host = host or cf["host"] or "ftp-private.ncbi.nlm.nih.gov"
        self.ncbiUser = cf["ncbiUser"]
        self.ncbiPass = cf["ncbiPass"]
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
        self.sars_cov_2_diag_pcr_ct_values = cf.get("sars_cov_2_diag_pcr_ct_values",{})
        self.skip_genbank = cf["skip_genbank"]
        self.sra_presets = cf["sra"]
        self.allowed_schemes = cf["allowed_schemes"]
        self.protocol_scheme = cf["protocol_scheme"]
        self.scheme = cf["scheme"]
        self.submitted_samples_file = cf.get("submitted_samples")
        self.ftp = None
        self.test_dir = test_dir if test_dir==True else cf.get("test_dir",test_dir)
        self.test_mode = test_mode if test_mode==True else cf.get("test_mode",test_mode)

        # set plate-specific details
        self.fastq_dir = Path(fastq_dir)
        self.seq_report = Path(seq_report)
        self.plate = plate
        self.outdir = Path(outdir)
        self.exclude_file = self.outdir / "samples2exclude.txt"
        self.barcode_map = None if barcode_map==None else Path(barcode_map) # only required for file prep
        self.report_dir = self.outdir / "bs_sra_reports"
        self.gisaid_log = gisaid_log
        self.primer_map = primer_map
        self.primer_scheme = primer_scheme
        self.fasta = fasta

        self.genbank = dict(
            name = "genbank",
            tsv = outdir / "genbank-metadata.tsv",
            xml_file = outdir / "genbank.xml",
            xml_text = None,
            df = None,
        )
        self.biosample = dict(
            name = "biosample",
            tsv = outdir / "biosample_attributes.tsv",
            xml_text = None,
            xml_file = outdir / "sra_biosample.xml",
            df = None,
        )
        self.sra = dict(
            name = "sra",
            tsv = outdir / "sra-metadata.tsv",
            xml_text = None,
            xml_file = outdir / "sra_biosample.xml",
            df = None,
        )
        self.merged = dict(
            tsv = outdir / "merged-metadata.tsv",
            df = None,
        )

        # TODO: rename spreadsheet as seq_report_df?
        self.spreadsheet = self._get_seq_report_df()
        self.barcode_df = self._prep_barcode_df()
        self.gisaid = self._prep_gisaid_df()
        self.gisaid_submitted_samples = None

    ### --- vvv --- miscellaneous --- vvv --- ###

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
        if self.seq_report == None: raise Exception("`seq_report` must be provided")
        spreadsheet = pd.read_csv(self.seq_report).rename(columns={"primer_scheme":"Primer Scheme","primer scheme":"Primer Scheme",})
        # if using biosample with headers beginning with underscores
        if hasUnderscoredHeaders(spreadsheet.columns):
            new_cols = {col:col.lstrip("_") for col in spreadsheet.columns if col.startswith("_")}
            spreadsheet = spreadsheet.rename(columns=new_cols)
            spreadsheet = spreadsheet.dropna(how='all')
        else:
            # if using initial EML rather than sequencing report:
            if "Seq ID" in spreadsheet.columns:
                spreadsheet = spreadsheet[[col for col in spreadsheet.columns if col != "Sample #"]]
                spreadsheet = spreadsheet.rename(columns={"Seq ID":"Sample #"})
            # weed out controls
            if "Sample #" in spreadsheet.columns:
                spreadsheet = spreadsheet[spreadsheet["Sample #"].str.contains(self.controls)==False]
            if "Test date" in spreadsheet.columns:
                spreadsheet["Test date"] = pd.DatetimeIndex(spreadsheet['Test date']).strftime('%Y-%m-%d')
            if 'Sample ID' in spreadsheet.columns:
                spreadsheet['Sample ID'] = spreadsheet['Sample ID'].astype(str)
        # spreadsheet['isolate'] = "SARS-CoV-2/human/USA/" + spreadsheet["Sample ID"].astype(str) + "/" + pd.DatetimeIndex(spreadsheet['Test date']).strftime('%Y')
        return spreadsheet

    def _prep_barcode_df(self):
        """Creates df from `barcode_map`"""

        if self.barcode_map == None:
            return pd.DataFrame()
        bc = pd.read_csv(self.barcode_map,sep="\t",names=["barcode","sample_name"])
        if self.controls and type(self.controls)==str:
            bc = bc[bc["sample_name"].str.contains(self.controls)==False]
        return bc

    def _prep_gisaid_df(self):
        """Creates df from `gisaid_log` - maps Virus name to accession"""

        # merge in gisaid accession details & update sample_name to include our location-based prefix
        if self.gisaid_log == None:
            return pd.DataFrame()
        else:
            gisaid = pd.read_csv(self.gisaid_log, sep=';', header=None, names=["gisaid_virus_name", "gisaid_accession"])
            for col in gisaid.columns: # remove extra spaces in values
                gisaid[col] = gisaid[col].apply(lambda x: str(x).strip())
                gisaid["sample_name"] = gisaid["gisaid_virus_name"].apply(lambda x: x.split("/")[2].replace(f"{self.affiliation['sub']}-",""))
            return gisaid[~gisaid["gisaid_virus_name"].str.contains("submissions")]

    def _ensure_new_names_only(self,df):
        """Raises Exception if any 'sample_name' in `df` was already-submitted"""

        if not self.submitted_samples_file: return
        if not Path(self.submitted_samples_file).exists():
            print(f"Couldn't find file: {self.submitted_samples_file}.\nNot checking for previously submitted samples.")
            return
        submitted = pd.read_csv(self.submitted_samples_file,header=None,names=["already_submitted"])["already_submitted"].unique()
        attempted_duplicate = df[df["sample_name"].isin(submitted)]
        if len(attempted_duplicate) > 0:
            raise Exception(f"You are attempting to submit a sample with a 'sample_name' that you have already submitted.\n\n{attempted_duplicate}\n\n * If you are sure you want to submit this sample, either comment out the variable `submitted_samples` in your `config_file` or remove that sample name from the file ({self.submitted_samples_file}). \n * If this sample should not be submitted, add it to {self.exclude_file}")

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

        if type(self.biosample["df"])==pd.DataFrame: return self.biosample["df"]
        # ^^ return existing df or vv create anew
        spreadsheet = self.spreadsheet
        # get variables from spreadsheet and set column headers to match those in template
        starting_cols = [x.strip() for x in self.biosample_presets.get("all_cols","").split(",") if x.strip() in spreadsheet.columns]
        if not starting_cols:
            starting_cols = [col for col in ["Sample #","Ct N gene","Test date","Sample ID","bioproject_accession","Primer Scheme"] if col in spreadsheet.columns]
            biosample_df = spreadsheet[starting_cols].rename(columns={"Test date":"collection_date","Sample #":"sample_name","Sample ID":"sample_title","Ct N gene":"sars_cov_2_diag_pcr_ct_value_1"})
        else:
            biosample_df = spreadsheet

        # add derived attributes
        # ensure these are added to the df and list of expected columns
        ct_cols = []
        for i,col in self.sars_cov_2_diag_pcr_ct_values.items():
            biosample_df[f'sars_cov_2_diag_pcr_ct_value_{i}'] = spreadsheet[col]
            ct_cols.extend([f'sars_cov_2_diag_pcr_ct_value_{i}',f'sars_cov_2_diag_gene_name_{i}'])

        # merge in gisaid accessions
        # biosample["sample_name_short"] = biosample['sample_name'].astype(str).apply(lambda x: x.split("-")[-1])
        # biosample = pd.merge(biosample,gisaid,on="sample_name_short",how='outer')
        if not self.gisaid.empty:
            # biosample = pd.merge(biosample,gisaid,on="sample_name_short",how='right')
            biosample_df = pd.merge(biosample_df,self.gisaid,on="sample_name",how='right')
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
        # print("reqs:",required_bs_attr)
        all_bs_cols = [x.strip() for x in self.biosample_presets.get("all_cols","").split(",") if x]
        # print("abc:",all_bs_cols)
        if not all_bs_cols:
            all_bs_cols = ['sample_name','sample_title','bioproject_accession','organism','collected_by','collection_date','geo_loc_name','host','host_disease','isolate','isolation_source','antiviral_treatment_agent','collection_device','collection_method','date_of_prior_antiviral_treat','date_of_prior_sars_cov_2_infection','date_of_sars_cov_2_vaccination','exposure_event','geo_loc_exposure','gisaid_accession','gisaid_virus_name','host_age','host_anatomical_material','host_anatomical_part','host_body_product','host_disease_outcome','host_health_state','host_recent_travel_loc','host_recent_travel_return_date','host_sex','host_specimen_voucher','host_subject_id','lat_lon','passage_method','passage_number','prior_sars_cov_2_antiviral_treat','prior_sars_cov_2_infection','prior_sars_cov_2_vaccination','purpose_of_sampling','purpose_of_sequencing']
            all_bs_cols.extend(ct_cols)
            all_bs_cols.extend(['sequenced_by','vaccine_received','virus_isolate_of_prior_infection','description'])
        if ignore_dates:
            required_bs_attr = remove_item("collection_date",required_bs_attr)
            all_bs_cols = remove_item("collection_date",all_bs_cols)
        # print("spreadsheet cols:\n",spreadsheet.columns)
        # bs_temp_cols = ["sample_name_short",]
        # actual_cols = [col for col in all_bs_cols+bs_temp_cols if col in biosample.columns]
        # extra_cols = [col for col in biosample_df.keys() if col not in set(all_bs_cols)]
        # print("extras:",extra_cols)
        actual_cols = [col for col in all_bs_cols if col in biosample_df.columns]
        # print(actual_cols)
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
                print("\nWARNING:\nSamples exist in barcode file that aren't found in BioSample file:\n\t" + str(extras))
                print(f"Any sample can be skipped by writing it in a line by itself in a file called:\n"
                    f"{self.exclude_file}\n"
                    f"If {'these samples' if len(extras)!=1 else 'this sample'} should be excluded from submissions, you can use this command:")
                for sample in extras:
                    print(f"echo '{sample}' >> '{self.exclude_file}'")
                warn("")
        # ensure all samples present have not been previously submitted
        self._ensure_new_names_only(biosample_df)
        print("BioSample")
        # print(biosample[["sample_name_short","bioproject_accession"]])
        print(biosample_df)
        self.biosample["final_cols"] = actual_cols

        # print(biosample_df.columns)
        self.biosample["df"] = biosample_df
        return biosample_df

    # def seekAccessionsDict(self):
    def get_accessions(self,as_dict=False,as_df=False,biosample_accessions:Path=None):
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
            # verify reports exist
            if not Path(self.report_dir).exists():
                raise FileNotFoundError(f"BioSample report directory not found:\n{self.report_dir}")
            accessions = {}
            # search for report with accessions
            for file in sorted([f for f in self.report_dir.glob("*/report.*") if f.name != "report.xml"],reverse=True):
                # print(f"checking {file} for accessions")
                report = Report(file)
                if report.biosamplesOk():
                    accessions_in_file = report.getAccessionDict(by_sample_name=True)
                    if accessions_in_file:
                        print("found accessions in",file)
                        # print(accessions_in_file)
                        accessions.update(accessions_in_file)
            if as_df:
                return pd.DataFrame(accessions.items(),columns=["Sequence_ID","BioSample"])
            elif as_dict:
                return accessions

    def add_biosamples_2_genbank(self,biosample_accessions=None):
        """Adds BioSample accessions to GenBank TSV"""

        # read in dfs
        print("reading in genbank")
        genbank_df = pd.read_csv(self.genbank["tsv"],sep="\t")
        g_cols = genbank_df.columns
        genbank_df = genbank_df.drop(columns="BioSample")
        print("reading in accessions")
        accessions_df = self.get_accessions(as_df=True,biosample_accessions=biosample_accessions)

        print("merging")
        genbank_df = genbank_df.merge(accessions_df,on="Sequence_ID",how="outer")[g_cols]
        genbank_df = genbank_df[genbank_df['note'].notna()]
        print(genbank_df)

        print("overwrite/save old genbank df to self")
        self.genbank["df"] = genbank_df

    def get_gisaid_submitted_samples(self):
        """Returns list of samples that were submitted to gisaid (based on `gisaid_log`) or "all" """

        if self.gisaid_submitted_samples == None:
            if not self.gisaid.empty:
                return list(self.gisaid["sample_name"])
            else: return "all"
        else: return self.gisaid_submitted_samples

    def _prep_genbank_df(self,add_biosample=False,biosample_accessions=None,ignore_dates=False):
        """Creates GenBank df or adds BioSample accessions to df from existing GenBank TSV

        Args:
            add_biosample (bool, optional): A flag to locate and add BioSample accessions to existing GenBank TSV. Defaults to False.
            biosample_accessions (str | Path, optional): path to a file downloaded from NCBI containing BioSample accessions
                if not provided, accessions will be sought from report.xml files
            ignore_dates (bool, optional): A flag to drop date columns from df
        """

        if add_biosample:
            self.add_biosamples_2_genbank(biosample_accessions)
            return
        
        cols = ['sample_name','organism','geo_loc_name', 'host', 'isolate', 'collection_date', 'isolation_source', 'bioproject_accession', 'gisaid_accession']
        if ignore_dates: cols = remove_item('collection_date',cols)
        genbank_df: pd.DataFrame = self.biosample["df"].copy()[cols]
        genbank_df['biosample_accession'] = 'Missing'
        genbank_df['gisaid_accession'] = genbank_df['gisaid_accession'].apply(lambda x: 'GISIAD accession: ' + str(x).strip() if str(x)!='missing' else "")
        submittable = self.get_gisaid_submitted_samples()
        # filter to gisaid-only
        if not submittable=="all":
            genbank_df = genbank_df[genbank_df["sample_name"].isin(submittable)]
        new_cols = cols + ['biosample_accession']
        genbank_df = genbank_df[new_cols]
        genbank_df = genbank_df.rename(columns={'sample_name':'Sequence_ID','geo_loc_name':'country', 'host':'host', 'isolate':'isolate', 'collection_date':'collection-date', 'isolation_source':'isolation-source', 'biosample_accession':'BioSample', 'bioproject_accession':'BioProject', 'gisaid_accession':'note'})
        print("GENBANK")
        print(genbank_df)
        # print(genbank.columns)
        self.genbank["df"] = genbank_df
        self.genbank["final_cols"] = list(genbank_df.columns)
        return genbank_df
    
    def _offer_skip_option(self,sample_name):
        return textwrap.dedent(f"""If submission of this sample should be skipped, write it in a line by itself in a file called:
            {self.exclude_file}\n"
            You can use this command:\n"
            echo "{sample_name}" >> {self.exclude_file}'""")

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
            FileNotFoundError(f"Too many potential fastq files for sample '{sample_name}' in {self.fastq_dir}")
        else:
            FileNotFoundError(f"Can't find expected fastq file for {sample_name} in {self.fastq_dir}\n{self._offer_skip_option(sample_name)}")
        if file_path.exists():
            return file_path.name
        else:
            raise FileNotFoundError(f"expected fastq file {file_path} does not exist.\n{self._offer_skip_option(sample_name)}")
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
            else:
                raise AttributeError(f"'Primer Scheme' must be a field in the `primer_map` file '{self.primer_map}'. Alternatively, provide `primer_scheme` as an argument.")
        # if no primer_map or 'Primer Scheme' column isn't found therein, use default for all samples
        elif self.primer_scheme:
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

        # get repeated fields from BioSample
        sra_df = self.biosample["df"].copy()
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
        print("SRA")
        print(sra_df)
        self.sra["df"] = sra_df
        self.sra["final_cols"] = actual_cols
        return sra_df

    # def merge_dfs(self,data,df_names):
    def _prep_merged_df(self):
        """Creates combined df of all samples/attributes"""

        biosample_df, sra_df, genbank_df = [getattr(self,dataset)["df"] for dataset in ("biosample","sra","genbank")]
        biosample_df:pd.DataFrame
        sra_df:pd.DataFrame
        genbank_df:pd.DataFrame
        # create merged df (before genbank['sample_title'] gets dropped)
        s1=set(biosample_df.columns)
        s2=set(sra_df.columns)
        # sra_extras = [col for col in s1.intersection(s2) if col != 'sample_name_short']
        # merged = biosample.merge(sra.drop(columns=sra_extras),on='sample_name_short',how='outer')
        # sra_extras = [col for col in s1.intersection(s2) if col != "sample_name"]
        sra_extras = remove_item('sample_name',s1.intersection(s2))
        # print(biosample.columns)
        # print(sra.drop(columns=sra_extras).columns)
        merged_df = biosample_df.merge(sra_df.drop(columns=sra_extras),on='sample_name',how='outer')
        # print(merged_df['sample_name'])
        # exit()
        if not genbank_df.empty:
            s3=set(genbank_df.columns)
            # genbank_extras = [col for col in s1.intersection(s3) if col != 'sample_name_short']
            # merged = merged.merge(genbank.drop(columns=genbank_extras),on='sample_name_short',how='outer')
            # genbank_extras = list(s1.intersection(s3))
            # genbank_extras = [col for col in s1.intersection(s3) if col != 'sample_name']
            genbank_extras = remove_item('sample_name',s1.intersection(s3))
            # print(genbank_df['sample_name'])
            # print(genbank_extras)
            # exit()
            merged_df = merged_df.merge(genbank_df.drop(columns=genbank_extras),left_on='sample_name',right_on='Sequence_ID',how='outer')
        # print(merged)
        self.merged["df"] = merged_df
        return merged_df

    def prep_dfs(self,add_biosample=False,ignore_dates=False,biosample_accessions=None):
        """Prepares spreadsheets for NCBI submission portal (https://submit.ncbi.nlm.nih.gov/)

        Args:
            sra (bool, optional): A flag to create SRA TSV. Defaults to True.
            biosample (bool, optional): A flag to create BioSample TSV. Defaults to True.
            genbank (bool, optional): A flag to create GenBank TSV. Defaults to True.
            add_biosample (bool, optional): A flag to add biosample accessions to GenBank TSV. Defaults to False.
            ignore_dates (bool, optional): A flag to drop date columns. Defaults to False.
        """

        if add_biosample:
            self._prep_genbank_df(add_biosample)
        else:
            self._prep_biosample_df(ignore_dates)
            self._prep_sra_df()
            self._prep_genbank_df(biosample_accessions=biosample_accessions,ignore_dates=ignore_dates)
            self._prep_merged_df()

    ###############  ^^^^  tsv_prep  ^^^^  ###############
    ###############  vvvv  xml_prep  vvvv  ###############

    ###############  ^^^^  xml_prep  ^^^^  ###############
    ###############  vvvv  zip_prep  vvvv  ###############

    def generate_updated_fasta_records(fasta):
        "returns generator of fasta records with only first part of header (drops /medaka something-or-other)"
        for record in SeqIO.parse(fasta,"fasta"):
            # update id
            record.id = record.id.split("/")[0]
            yield record

    def get_sample_name(virus_name):
        "converts GISAID Virus name to our typical sample identifier"
        return virus_name.split("/")[2].replace("NC-","") # assumes state is NC # TODO: generalize

    def getAllowedSamples(data,gisaid_log):
        "uses only samples in gisaid_log (if provided) or else previously made genbank TSV as allowed samples from provided fasta"
        # acquire collection of allowed samples from gisaid logfile, if provided
        # gisaid_log = getattr(args,"gisaid_log",None)
        if gisaid_log != None:
            # allowed_samples = set(pd.read_csv(gisaid_log,header=None,names=["gisaid_virus_name","Accession ID"],sep=";").dropna()["gisaid_virus_name"].apply(get_sample_name))
            allowed_samples = set(getGisaidDF(gisaid_log)["gisaid_virus_name"].apply(get_sample_name))
        else:
            # assume original TSV from file_prep step is filtered to gisaid-submitted (or otherwise desired) samples
            allowed_samples = set(data["GenBank"].df["Sequence_ID"].unique())
        return allowed_samples

    def write_genbank_fasta(outHandle,data,fasta):
        "writes GenBank fasta to provided outHandle"
        allowed_samples = getAllowedSamples(data,gisaid_log)
        # print(allowed_samples)
        for record in generate_updated_fasta_records(fasta):
            if record.id in allowed_samples:
                outHandle.write(f">{record.id}\n")
                # trim leading/trailing Ns
                record.seq = record.seq.strip("N")
                outHandle.write(f"{record.seq}\n")

    def ensureCorrectBioProject(tsv_str,test_dir):
        "if this is a test submission, replace all BioProject accessions with one that exists in the test submission database"
        if not test_dir: return tsv_str
        row_list = []
        for row in tsv_str.split("\n"):
            listed_row = []    
            for col in row.split("\t"):
                if "PRJNA" in col:
                    col = "PRJNA553747"
                listed_row.append(col)
            row_list.append("\t".join(listed_row))
        return "\n".join(row_list)

    def write_genbank_zip(self,data,outdir:Path,test_dir,template,plate):
        "writes zipfile for GenBank submission"
        # documentation and sample files are here: https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/genbank/SARS-CoV-2/
        with ZipFile(outdir.joinpath("genbank.zip"), "w") as zfh:
            # writing from stringio avoids any py2/3 unicode trouble (-Multisub)
            print("\twriting fasta")
            with StringIO() as memFh:
                self.write_genbank_fasta(memFh,data,fasta)
                zfh.writestr("seqs.fsa", memFh.getvalue())

            print("\twriting tsv --> src")
            with outdir.joinpath("genbank-metadata.tsv").open() as tsv:
                tsv_string = self.ensureCorrectBioProject(tsv.read(),test_dir)
                zfh.writestr("seqs.src", tsv_string)
            # with StringIO() as memFh:
            #     writeMetaTsv(sourceTags, memFh)
            #     zfh.writestr("seqs.src", memFh.getvalue())

            print("\twriting template")
            # template = getattr(args,"template",None) # priority goes to template path from args
            if template==None:
                template = cf.get("template",None) # seek template path from `config_file`
            if not template:
                # TODO: generate template from `config_file`
                raise Exception("Path to `template` file must be provided in `config_file` or args.")
            with Path(template).open() as sbt, StringIO() as memFh:
                for line in sbt:
                    if "Submission Title:None" in line:
                        line = line.replace(":None",f":{plate} GenBank")
                    memFh.write(line)
                zfh.writestr("seqs.sbt",memFh.getvalue())
                # zfh.writestr("seqs.sbt",sbt.read())
            # zfh.writestr("seqs.sbt", makeTemplate())

    ###############  ^^^^  zip_prep  ^^^^  ###############

    def write_genbank_xml(self):
        """Writes out XML for GenBank submission"""

        xml_maker = GenBank_Submission(self)
        xml_maker.write_xml("genbank.xml")

    def write_sra_biosample_xml(self):
        """Writes out XML for SRA and/or BioSample submission"""

        xml_maker = SRA_BioSample_Submission(self)
        xml_maker.write_xml("sra_biosample.xml")

    def write_presubmission_metadata(self,sra_only=False,keep_tsvs=True,keep_xmls=True):
        """Prepares TSV and XML files

        This should be the go-to method for preparing and writing out files for initial submission to BioSample and/or SRA

        Args:
            action (str: "file_prep" | "add_biosamples"): 
                "file_prep": Creates initial TSV and XML files
                "add_biosamples": adds biosample accessions from previous submission to genbank tsv
            keep_tsvs (bool, optional): If False, TSVs will be deleted. Defaults to True.
            keep_xmls (bool, optional): If False, XMLs will be deleted. Defaults to True.
        """

        ignore_dates = True if sra_only else False

        # prepare/write TSVs
        self.prep_dfs(ignore_dates=ignore_dates)
        # final_dfs = {} # NOTE: not necessary?
        for dataset in (self.sra, self.biosample, self.genbank):
            # final_df = finalize_df(df=dataset["df"],final_cols=dataset.get("final_cols"))
            # # final_dfs[dataset["name"]] = final_df # NOTE: not necessary?
            # df_2_tsv(df=final_df,outfile=dataset["tsv"],name=dataset["name"])
            # exit()
            df_2_tsv(
                df=finalize_df(df=dataset["df"],final_cols=dataset["final_cols"]),
                outfile=dataset["tsv"],name=dataset["name"]
            )

        # prepare/write XML
        print("writing sra/biosample xml")
        self.write_sra_biosample_xml()
        self.write_genbank_xml()

        # # remove unwanted text files
        # if not keep_tsvs:
        #     for dataset in (self.sra, self.biosample, self.genbank):
        #         tsv:Path = dataset["tsv"]
        #         tsv.unlink()

        # if not keep_xmls:
        #     for dataset in (self.sra, self.biosample, self.genbank):
        #         xml:Path = dataset["xml_file"]
        #         xml.unlink()




    ### --- ^^^ --- file prep --- ^^^ --- ###

    ### --- vvv --- ncbi interaction --- vvv --- ###
    # NOTE: TODO: add in replaceTestAccessions() before submitting