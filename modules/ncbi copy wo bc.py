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
        self.test_dir = cf.get("test_dir",test_dir)
        self.test_mode = cf.get("test_mode",test_mode)

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

        # TODO: rename spreasdheet as seq_report_df?
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

    ### --- ^^^ --- miscellaneous --- ^^^ --- ###
    ### --- vvv --- file prep --- vvv --- ###

    def _get_seq_report_df(self):
        """Creates df from `seq_report`"""

        # local sequencing report data
        if self.seq_report == None: raise Exception("`seq_report` must be provided")
        spreadsheet = pd.read_csv(self.seq_report)
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
            raise Exception(f"You are attempting to submit a sample with a 'sample_name' that you have already submitted.\n\n{attempted_duplicate}\n\n * If you are sure you want to submit this sample, either comment out the variable `submitted_samples` in your `config` file or remove that sample name from the file ({self.submitted_samples_file}). \n * If this sample should not be submitted, add it to {self.exclude_file}")

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

        print(df["bioproject_accession"])
        default_production_accession = self.bioproject_presets.get("bioproject_accession")   # vvv a random test-dir bioproject accession
        default_test_accession = self.bioproject_presets.get("bioproject_test_accession","PRJNA553747")
        if self.test_dir:
            if not isBioProjectAccession(default_test_accession):
                raise ValueError(f"Default test BioProject accession '{default_test_accession}' in config '{self.config_file}' does not fit the standard format (PRJ[D|E|N]xxxxxx or PSUBxxxxxx).")
            df["bioproject_accession"] = default_test_accession
        else:
            if not isBioProjectAccession(default_production_accession):
                raise ValueError(f"Default BioProject accession '{default_production_accession}' in config '{self.config_file}' does not fit the standard format (PRJ[D|E|N]xxxxxx or PSUBxxxxxx).")
            if not "bioproject_accession" in df.columns:
                df["bioproject_accession"] = default_production_accession
            else:
                df = df.fillna(default_production_accession)
                # ensure all accessions look like accessions
                checkBioProjectAccessions(df["bioproject_accession"])
        print(df["bioproject_accession"])
        exit()
        return df

    def _prep_biosample_df(self,ignore_dates=False):
        """Creates BioSample df based on data in `seq_report` and `config`"""

        if type(self.biosample["df"])==pd.DataFrame: return self.biosample["df"]
        # ^^ return existing df or vv create anew
        spreadsheet = self.spreadsheet
        # get variables from spreadsheet and set column headers to match those in template
        starting_cols = [x.strip() for x in self.biosample_presets.get("all_cols","").split(",") if x]
        if not starting_cols:
            starting_cols = [col for col in ["Sample #","Ct N gene","Test date","Sample ID","bioproject_accession"] if col in spreadsheet.columns]
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
        bc_defaults = self.biosample_presets
        for k,col in bc_defaults.items():
            if not k in ("bioproject_accession","bioproject_test_accession"):
                biosample_df[k] = col

        # use test bioproject for all samples (if --test_dir==True)
        # use default bioproject for any samples that don't have one (if --test_dir==False)
        biosample_df = self._add_bioproject_if_needed(biosample_df)
        exit()

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
        extra_cols = [col for col in biosample_df.keys() if col not in set(all_bs_cols)]
        actual_cols = [col for col in all_bs_cols if col in biosample_df.columns]
        print(actual_cols)
        check_missing(actual_cols,required_cols=required_bs_attr,name="BioSample")
        biosample_df = biosample_df[actual_cols]
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
        self.biosample["df"] = biosample_df
        return biosample_df

    # def seekAccessionsDict(self):
    def getAccessions(self,as_dict=False,as_df=False,biosample_accessions:Path=None):
        """Returns dict or df mapping sequence id to biosample accession

        Args:
            as_dict (bool, optional): _description_. Defaults to False.
            as_df (bool, optional): _description_. Defaults to False.
            biosample_accessions (Path, optional): path to accession file. Defaults to None.
                If provided, uses file to determine accessions
                If not, recursively searches all report*.xml files in report_dir for one containing BioSample accessions

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
        accessions_df = self.getAccessions(as_df=True,biosample_accessions=biosample_accessions)

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
        genbank_df.rename(columns={'sample_name':'Sequence_ID','geo_loc_name':'country', 'host':'host', 'isolate':'isolate', 'collection_date':'collection-date', 'isolation_source':'isolation-source', 'biosample_accession':'BioSample', 'bioproject_accession':'BioProject', 'gisaid_accession':'note'}, inplace=True)
        print("GENBANK")
        print(genbank_df)
        # print(genbank.columns)
        self.genbank["df"] = genbank_df
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

        # file_path:Path = self.fastq_dir / f"{sample_name}_{row['barcode']}.fastq" #[0]
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
        """Adds primer scheme column to SRA df if possible/needed"""
        if "Primer Scheme" in df.columns:
            return df
        # attempt to get scheme from primer map file
        if self.primer_map:
            primers = pd.read_csv(self.primer_map).rename(columns={"Sample #":"sample_name","Seq ID":"sample_name"})
            if "Primer Scheme" in primers.columns:
                primers = primers[["sample_name","Primer Scheme"]]
                print(primers)
                df = df.merge(primers,on="sample_name")
                # print(df.columns)
                # print(df)
                # print("Must check this code (in `add_primer_schemes()`) before submitting anything with primer scheme col in primer_map")
                # exit(1)
                return df
        # if no primer_map or scheme column not found therein, use default for all samples
        if self.primer_scheme:
            df["Primer Scheme"] = self.primer_scheme
            return df
        else: raise Exception("A 'Primer Scheme' field must be in your seq_report file or else \n one of the arguments `--primer_map` or `--primer_scheme` must be provided")

    def _addFilenames(self,sra_df:pd.DataFrame):
        """Adds paths to filenames for each sample (pair reads get extra filename column)"""

        if "illumina" in self.sra_presets["platform"].lower():
            sra_df[["filename","filename2"]] = sra_df['sample_name'].apply(self._get_paired_fastq_files,result_type="expand")
        else:
            sra_df["filename"] = sra_df["sample_name"].apply(self._get_fastq_file)
        return sra_df

    def _prep_sra_df(self):
        """Returns DataFrame of SRA attributes derived from BioSample data or retrieved from `config`"""

        # get repeated fields from BioSample
        sra_df = self.biosample["df"][['bioproject_accession','sample_name','sample_title']].copy()
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

        # add `config` defaults/overrides
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
        sra_df = sra_df[actual_cols]
        check_missing(actual_cols,required_sra_cols,"SRA")
        print("SRA")
        print(sra_df)
        self.sra["df"] = sra_df
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

    def get_title(self):
        """Returns unique title line for plate"""

        if self.plate: return f"""\
        <Title>{self.plate} - {self.biosample_presets["isolation_source"].lower()} - sars-cov-2</Title>"""
        else: return ""

    def ncbi_header(self,data_type,test_dir:bool):
        """Returns lab-specific header for any NCBI xml file
        
        Args:
        data_type (str): name of submission type ("BS+SRA","GB",...)
        test_dir (bool): A flag to prepare a test (non-production) submission
        """

        ## TODO: allow date specifications in `config`?
        # next_week = datetime.datetime.now() + datetime.timedelta(days=7)
        # release = next_week.strftime('%Y-%m-%d')
        today = datetime.datetime.now()
        release = today.strftime('%Y-%m-%d') #TODO: add this to `config`/parameters
        test = "test " if test_dir==True else ""
        return f"""<?xml version="1.0"?>
    <Submission>
    <Description>
        {self.get_title().lstrip()}
        <Comment>SARS-CoV-2 {test}submission - {self.plate} {data_type}</Comment>
        <Organization type="center" role="owner">
        <Name>{self.affiliation['div']}</Name>
        <Contact email="{self.email}">
            <Name>
            <First>{self.contact[1]}</First>
            <Last>{self.contact[0]}</Last>
        </Name>
        </Contact>
        </Organization>
        <Hold release_date="{release}"/>
    </Description>"""

    def add_BioSample_xml(self,row,biosample_spuid,cf,plate,test_dir): #TODO: edit all this
        "returns xml action to submit given sample to BioSample"
        # package options and details: https://www.ncbi.nlm.nih.gov/biosample/docs/packages/
        # SARS-CoV-2.cl.1.0 is assumed, currently but this could be moved to config
        package = cf["biosample"]["package"]
        xml_text_list = []
        xml_text_list.append(f"""\
    <Action>
        <AddData target_db="BioSample">
            <Data content_type="XML">
            <XmlContent>
                <BioSample schema_version="2.0">
                <SampleId>
                    {biosample_spuid}
                </SampleId>
                <Descriptor>
                    <Title>{plate} BioSample</Title>
                </Descriptor>
                <Organism>
                    <OrganismName>{cf["biosample"]["organism"]}</OrganismName>
                </Organism>""")
        if not cf['bioproject']['create_new']: xml_text_list.append(f"""\
                <BioProject>
                    <PrimaryId db="BioProject">{row["bioproject_accession"]}</PrimaryId>
                </BioProject>""")
        xml_text_list.append(f"""\
                <Package>{package}</Package>
                <Attributes>
                    {add_attributes(row,db="biosample",blankspace="                ").lstrip()}
                </Attributes>
                </BioSample>
            </XmlContent>
            </Data>""")
        if cf['bioproject']['create_new']: xml_text_list.append(f"""\
        <AttributeRefId name="BioProject">
            <RefId>
            {get_bioproject_spuid(cf)}
            </RefId>
        </AttributeRefId>""")
        xml_text_list.append(f"""\
        <Identifier>
            {biosample_spuid}
        </Identifier>
        </AddData>
    </Action>""")
        return "\n".join(xml_text_list)

    def add_attributes(row,db,attrs=[],blankspace=""):
        "returns xml text for adding an attribute to the requested db"
        if not attrs: attrs = [x for x in row.index if not "bioproject_accession" in x]
        attr_text_list = []
        for attr_name in attrs:
            attr = row[attr_name]
            # skip empty attributes
            if str(attr) == "nan": continue
            elif not str(attr).strip(): continue
            # add attributes with different spacing for different submissions
            if db == "sra":
                name_of_attr = "name"
                # these attributes are excluded (since they go elsewhere)
                if "filename" in attr_name or "filetype" in attr_name:
                    continue
            elif db == "biosample":
                name_of_attr = "attribute_name"
            else: warn(f"db {db} not accounted for")
            attr_text_list.append(f'{blankspace}<Attribute {name_of_attr}="{attr_name}">{attr}</Attribute>')
        return "\n".join(attr_text_list)

    def add_file(file,datatype="generic-data"):
        "returns xml text for adding a fastq file to an SRA sample"
        if str(file) == "nan":
            raise Exception("Filename looks like 'nan' - verify that all samples are labeled correctly in this plate's gisaid_uploader.log file")
        return f"""\
        <File file_path="{file}">
            <DataType>{datatype}</DataType>
        </File>"""

    def add_files_from_row(row):
        "returns xml text for adding all fastq files to an SRA sample"
        file_list = []
        for col in row.index:
            if type(col)==str and "filename" in col:
                file_list.append(add_file(row[col]))
        return "\n".join(file_list)

    def getBioBrojectLink(accession=None):
        "returns link to bioproject"
        if accession:
            return f'<PrimaryId db="BioProject">{accession}</PrimaryId>'
        else: # need to create BioProject #TODO: make sure this link is correct for BioProject creation
            print("getBioBrojectLink(None) - no accession given - code not ready for this")
            exit()
            # return f'<PrimaryId db="BioProject">{accession}</PrimaryId>'

    def link_bioProject(accession,cf): # this could be replaced by a function that only links biosample if not creating it as a seperate Action step
        "returns xml text for refencing a BioProject accession to a sample's BioSample submission"
        return f"""\
        <AttributeRefId name="BioProject">
            <RefId>
            {getBioBrojectLink(accession)}
            </RefId>
        </AttributeRefId>"""

    def link_bioproject_by_spuid(bioproject_spuid):
        "returns xml text referencing corresponding BioProject (for SRA submission)"
        return f"""\
        <AttributeRefId name="BioProject">
            <RefId>
            {bioproject_spuid}
            </RefId>
        </AttributeRefId>"""

    def link_bioSample_by_spuid(biosample_spuid):
        "returns xml text referencing corresponding BioSample (for SRA submission)"
        return f"""\
        <AttributeRefId name="BioSample">
            <RefId>
            {biosample_spuid}
            </RefId>
        </AttributeRefId>"""

    def add_sra_xml(row,sra_spuid,biosample_spuid,bioproject_accession,cf):
        "returns xml action to submit given sample to SRA"
        xml_text_list = []
        xml_text_list.append(f"""\
    <Action>
        <AddFiles target_db="SRA">
        {add_files_from_row(row).lstrip()}
        {add_attributes(row,db="sra",blankspace="      ").lstrip()}""")
        if cf["bioproject"]["create_new"] == False: xml_text_list.append(f"""\
        {link_bioProject(bioproject_accession,cf).lstrip()}""")
        xml_text_list.append(f"""\
        {link_bioSample_by_spuid(biosample_spuid).lstrip()}
        <Identifier>
            {sra_spuid}
        </Identifier>
        </AddFiles>
    </Action>""")
        return "\n".join(xml_text_list)

    def getBioProjectSpuid(cf):
        return f'<SPUID spuid_namespace="{cf["centerAbbr"]}">{cf["bioproject"]["spuid"]}</SPUID>'

    def add_BioProject_xml(cf):
        "returns xml action to create a new BioProject"
        bioProjectSpuid = get_bioproject_spuid(cf)
        # see https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/bioproject/bioproject.xsd?view=co
        return f"""\
    <Action>
        <AddData target_db="BioProject">
        <Data content_type="xml">
            <XmlContent>
            <Project schema_version="2.0">
                <ProjectID>
                {bioProjectSpuid}
                </ProjectID>
                <Descriptor>
                <Title>{cf["bioproject"]["spuid"]}</Title>
                <Description>
                    <p>{cf["bioproject"]["description"]}</p>
                </Description>
                <ExternalLink label="Website: {cf["bioproject"]["website_name"]}">
                    <URL>{cf["bioproject"]["url"]}</URL>
                </ExternalLink>
                <Relevance>
                    <Medical>Yes</Medical>
                </Relevance>
                </Descriptor>
                <ProjectType>
                <ProjectTypeSubmission sample_scope="{cf["bioproject"]["scope"]}">
                    <IntendedDataTypeSet>
                    <DataType>{cf["bioproject"]["dataType"]}</DataType>
                    </IntendedDataTypeSet>
                </ProjectTypeSubmission>
                </ProjectType>
            </Project>
            </XmlContent>
        </Data>
        <Identifier>
            {bioProjectSpuid}
        </Identifier>
        </AddData>
    </Action>"""

    def get_cols(db,outdir):
        "returns fields associated with the specified tsv (so the correct attributes can be pulled from the merged df)"
        dbs = {"sra":"sra-metadata","genbank":"genbank-metadata","biosample":"biosample_attributes"}
        return pd.read_csv(outdir / f"{dbs[db]}.tsv", sep="\t").columns

    def seekAccessionsDict(outdir:Path):
        "recursively searches all report*.xml files in report_dir for one containing BioSample accessions, returning them, if found"
        # verify reports exist
        report_dir = outdir/f"bs_sra_reports"
        if not Path(report_dir).exists():
            raise FileNotFoundError(f"BioSample report directory not found:\n{report_dir}")
        accessions = {}
        # search for report with accessions
        for file in sorted([f for f in report_dir.glob("*/report.*") if f.name != "report.xml"],reverse=True):
            # print(f"checking {file} for accessions")
            report = Report(file)
            report.biosamplesOk()
            if report.biosamplesOk():
                accessions_in_file = report.getAccessionDict(by_sample_name=True)
                if accessions_in_file:
                    print("found accessions in",file)
                    # print(accessions_in_file)
                    accessions.update(accessions_in_file)
        return accessions

    def get_sample_data(df,cf,outdir:Path,test_dir,sra_only,biosample_accessions,plate):
        "returns xml actions for sample's submission to BioSample and SRA"
        biosample_cols = get_cols("biosample",outdir)
        sra_cols = get_cols("sra",outdir)

        # determine accession to use
        if test_dir:
            bioproject_accession = cf["bioproject"]["bioproject_test_accession"].strip()
        else: bioproject_accession = cf["bioproject"]["bioproject_accession"].strip()
        if not "PRJNA" in bioproject_accession:
            yield add_BioProject_xml(cf)

        # sra_only = True #TEMP
        # sra_only = False #TEMP

        # print(args)
        # exit(1)
        if sra_only:
            # use provided accessions file
            if biosample_accessions:
                accession_dict = accessions_2_df(as_dict=True)
            # find accessions from submission report
            else:
                accession_dict = seekAccessionsDict(outdir)

        # generate xml text for each sample (either for biosample and sra or just sra)
        for i,row in df.iterrows():
            sra_link = f"""<SPUID spuid_namespace="{cf["centerAbbr"]}">{row["sample_name"]}_SRA</SPUID>"""
            if sra_only:
                biosample_link = f"""<PrimaryId db="BioSample">{accession_dict[row["sample_name"]]}</PrimaryId>"""
                xml_text = add_sra_xml(row[sra_cols],sra_link,biosample_link,bioproject_accession,cf) # temp
            else:
                biosample_link = f"""<SPUID spuid_namespace="{cf["centerAbbr"]}">{row["sample_name"]}_BioSample</SPUID>"""
                xml_text = add_BioSample_xml(row[biosample_cols],biosample_link,cf,plate,test_dir) + "\n" + add_sra_xml(row[sra_cols],sra_link,biosample_link,bioproject_accession,cf) # works
            yield xml_text
            # TODO: return spiuds for genbank?

    def write_sra_biosample_xml(cf,outdir:Path,test_dir,sra_only,biosample_accessions,plate):
        "writes out xml file for submitting all samples to BioSample and SRA"
        output = []
        output.append(self.ncbi_header("BS+SRA",plate,test_dir))
        df = pd.read_csv(outdir / "merged-metadata.tsv", sep="\t")
        for sample_data in get_sample_data(df,cf,outdir,test_dir,sra_only,biosample_accessions,plate):
            output.append(sample_data)
        # TODO: add references, as needed to genbank xml?
        output.append("</Submission>")
        with outdir.joinpath("sra_biosample.xml").open("w") as out:
            out.write("\n".join(output))

    def add_genbank_xml(cf,plate):
        "returns xml action for submitting all samples from plate to GenBank"
        return f"""\
    <Action>
        <AddFiles target_db="GenBank">
        <File file_path="genbank.zip">
            <DataType>genbank-submission-package</DataType>
        </File>
        <Attribute name="wizard">BankIt_SARSCoV2_api</Attribute>
        <Attribute name="auto_remove_failed_seqs">no</Attribute>
        <Identifier>
            <SPUID spuid_namespace="{cf["centerAbbr"]}">{plate}.sarscov2</SPUID>
        </Identifier>
        </AddFiles>
    </Action>"""

    def write_genbank_xml(cf,outdir:Path,plate,test_dir):
        "writes out xml file for submitting all samples to GenBank"
        output = "\n".join([
            ncbi_header(cf,"GB",plate,test_dir),
            add_genbank_xml(cf,plate),
            "</Submission>"
        ])
        with outdir.joinpath("genbank.xml").open("w") as out:
            out.write(output)

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

    def write_genbank_zip(cf,data,outdir:Path,test_dir,template,plate):
        "writes zipfile for GenBank submission"
        # documentation and sample files are here: https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/genbank/SARS-CoV-2/
        with ZipFile(outdir.joinpath("genbank.zip"), "w") as zfh:
            # writing from stringio avoids any py2/3 unicode trouble (-Multisub)
            print("\twriting fasta")
            with StringIO() as memFh:
                write_genbank_fasta(memFh,data,fasta)
                zfh.writestr("seqs.fsa", memFh.getvalue())

            print("\twriting tsv --> src")
            with outdir.joinpath("genbank-metadata.tsv").open() as tsv:
                tsv_string = ensureCorrectBioProject(tsv.read(),test_dir)
                zfh.writestr("seqs.src", tsv_string)
            # with StringIO() as memFh:
            #     writeMetaTsv(sourceTags, memFh)
            #     zfh.writestr("seqs.src", memFh.getvalue())

            print("\twriting template")
            # template = getattr(args,"template",None) # priority goes to template path from args
            if template==None:
                template = cf.get("template",None) # seek template path from `config`
            if not template: raise Exception("path to `template` file must be provided in `config` or args")
            with Path(template).open() as sbt, StringIO() as memFh:
                for line in sbt:
                    if "Submission Title:None" in line:
                        line = line.replace(":None",f":{plate} GenBank")
                    memFh.write(line)
                zfh.writestr("seqs.sbt",memFh.getvalue())
                # zfh.writestr("seqs.sbt",sbt.read())
            # zfh.writestr("seqs.sbt", makeTemplate())

    ###############  ^^^^  zip_prep  ^^^^  ###############

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
        for dataset in (self.sra, self.biosample, self.genbank):
            df_2_tsv(dataset["df"],dataset["tsv"],dataset['name'])

        # # prepare/write XMLs
        # self.prep_xmls()

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
