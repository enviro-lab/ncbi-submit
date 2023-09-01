# ncbi-submit
Submitting data to public databases is super important for publically funded laboratories, but it is not always a quick or intuitive process. `ncbi-submit` provides a simple and repeatable way to upload programmatic submissions to NCBI's SRA and GenBank with shared or unique BioProjects and BioSamples. Data can be uploaded as XML or zip files to either the Test or Production environments, and once there, the reports produced by NCBI can be analyzed to check on submission status and get BioSample accessions.

***
## Installation:
To install from PyPI in a virtual environment `.venv`:
```bash
python3 -m venv .venv
. .venv/bin/activate
pip install ncbi-submit
```
To install from conda (not yet set up) in a new environment `ncbi`:
```bash
conda create -n ncbi ncbi-submit
```

***
## Testing
Add NCBI credentials to file `./.login_credentials` or edit them in either:
* `./example/test.sh` or
* `./config/config.py`

To test creating all example files, run:
```bash
./example/test.sh
```
This script ^^^ could also be a good starting point for your own NCBI submission pipelines. Note: There are several blocks of code in there can be commented in/out, as needed.

***
## Usage

`ncbi_submit` is intended for use on the command line, but the class `ncbi.NCBI` can be imported and used within custom python scripts.

There are three main actions the script can do:
* `file_prep`: 
  * Prepares .tsv & .xml files for SRA, BioSample, & BioProject submissions
  * Used to prepare all files for initial submission to NCBI
  * To add in biosample accessions and prepare for GenBank submission, include the flag `prep_genbank`:
    * Prepares .zip, .sbt, & .tsv files for GenBank Submission
    * Used to add BioSample accessions from a BioSample submission for a GenBank submission
* `ftp` submission or checkup:
  * Interacts with NCBI's ftp host to do either of the following:
    * `submit` data to NCBI databases 
    * `check` on previous ftp submissions
    * `get-accessions` from all previous ftp submissions
* `example`:
  * Writes out example files for one or both of:
    * config.py file (tells `ncbi_submit` lots of important info)
    * template.sbt (used for genbank submission)

### Setup
The required parameters vary by which of the above actions you're attempting but at minimum require a `plate` and `outdir`. To limit the number of parameters required via command line, a `config` file must be used. When running from the command line, one of the three actions (`file_prep` or `ftp`) must be specified. With python, these are associated methods you may use on a single NCBI object.

#### Get example `config.py` file:
```bash
ncbi_submit example --config --outdir "nbci"
```

#### Python instantiation (not needed on command line):
Note: This is the minimum required info for preparing data. Other parameters may be necessary for more functionality or other tasks.  
```python
from ncbi_submit import ncbi_submit
ncbi = ncbi_submit.NCBI(
    fastq_dir = myFastqDir,
    seq_report = mySeqReport,
    plate = myPlate,
    outdir = myOutdir,
    config_file = myConfig,
    )
ncbi.write_presubmission_metadata()
```

### File Preperation
##### Shell:
```bash
ncbi_submit file_prep \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --seq_report "${SEQ_REPORT}" \
    --primer_map "${PRIMER_MAP}" \
    --primer_scheme "${SCHEME_VERSION}" \
    --outdir "${NCBI_DIR}" \
    --gisaid_log "${GENERIC_GISAID_LOG//PLATE/$PLATE}" \
    --fastq_dir ${FASTQS} \
    --plate "${PLATE}"
```
##### Python:
```python
ncbi.write_presubmission_metadata()
```

### File Submission
##### Shell:
```bash
ncbi_submit ftp submit \
    --db bs_sra \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --outdir "${NCBI_DIR}" \
    --fastq_dir "${FASTQS}"
# wait a while and try this to download reports and view submission status
ncbi_submit ftp \
    --check --db bs_sra \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --outdir "${NCBI_DIR}" 
```
##### Python:
```python
ncbi.submit(db="bs_sra")
# wait awhile and try this to download reports and view submission status
ncbi.check(db="bs_sra")
```

### GenBank submission
(NOTE: not fully tested)
To link your fasta in GenBank to the associated reads, you'll want to add in the BioSample accessions before submitting.
* Acquire BioSample accessions via one of these methods:
  * download accessions.tsv file from NCBI and then use `ncbi_submit`
    * (Do this if you submitted to BioSample via NCBI's Submission Portal)
  * use `ncbi_submit` for everything
    * (Do this to avoid manual uploads via NCBI's Submission Portal)
##### Shell:
```bash
# dowload report.xml files to get accesssions from
ncbi_submit ftp check \
    --db ${DB} \
    --outdir "${NCBI_DIR}" \
    --config "${NCBI_CONFIG}" \
    -u "${ncbi_username}" \
    -p "${ncbi_password}" \
    --plate "${PLATE}" \
    --fastq_dir "${FASTQS}"

# add accessions to genbank.tsv
ncbi_submit --prep_genbank \
    --outdir "${NCBI_DIR}" \
    --config ${NCBI_CONFIG} \
    --fasta "${GENERIC_CONSENSUS//PLATE/$PLATE}" \
    --plate "${PLATE}"

# submit to GenBank (NOTE: db='gb')
ncbi_submit ftp submit \
    --db gb \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --outdir "${NCBI_DIR}" \
    --fastq_dir "${FASTQS}"
```
##### Python:
```python
# dowload report.xml files to get accesssions
ncbi.check(db="bs_sra")
# prepare genbank submission files and submit
ncbi.submit(db="gb")

## or

# files can also be prepared without submitting via:
ncbi.write_genbank_submission_zip()
```

***
### Check Submission Status
Wait awhile (10+ minutes) for NCBI to start processing the submission. Then run this to download reports and view submission status.
This works for whichever db you want to check on. If not specified, you'll get results on all submitted dbs.

##### Shell:
```bash
# check GenBank submission status (NOTE: db='gb')
ncbi_submit ftp check \
    --db gb \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --outdir "${NCBI_DIR}"
```
##### Python:
```python
# check GenBank submission status (NOTE: db='gb')
ncbi.check(db='gb')
```

### Download all reports and get accessions
To acquire the accessions for all samples submitted via ftp under your group's account, `ncbi_submit` can download all xml report files and parse out the accession details. A directory will be created in `outdir` containing all submission-specific directories, each containing its report files. The `-f` or `--files` flag allows the use of a list of report files to parse. If provided, those files will be parsed for accession details rather than downloading the latest report files. The database can be specifed to indicate which accessions are desired and yield csvs (for the BioProject associated with your current `config` file) at `<outdir>/accessions_<bioproject>.csv` with the following fields:
| database | fields |
|-|-|
| 'bs_sra' | sample_name, BioSample, SRA |
| 'bs' | sample_name, BioSample |
| 'sra' | sample_name, SRA |

#### Get accessions by downloading report.xml files
##### Shell:
```bash
ncbi_submit ftp get-accessions \
    --db "bs_sra" \
    --config "${NCBI_CONFIG}" \
    --outdir "${REPORT_DIR}" \
    -u "${ncbi_username}" \
    -p "${ncbi_password}" \
```
##### Python:
```python
ncbi.get_all_accessions(db="bs_sra)
```

#### Get accessions from list of report.xml files
##### Shell:
```bash
ncbi_submit ftp get-accessions \
    --db "bs_sra" \
    --config "${NCBI_CONFIG}" \
    --outdir "${REPORT_DIR}" \
    -u "${ncbi_username}" \
    -p "${ncbi_password}" \
    -f s1/report.xml s2/report.xml
```
##### Python:
```python
ncbi.get_all_accessions(db="bs_sra",report_files=["file1", "file2"])
```


## Updating samples that have already been submitted
### Fastq read updates
If you want to update the reads for a sample you've already submitted, you must do the followind:
1. Email nlm-support@nlm.nih.gov and supply them with a list of SRA runs to suppress.
2. Once suppressed, you can upload a new version of the sample where the `submission.xml`
  * references the BioSample (rather than submitting a new BioSample block) and
  * has a new, unique SPUID for the SRA action block.

The `submission.xml` can be prepared as shown below and then submitted as discussed previously in [File Submission](#file-submission). Whereas normally an error would occur if a previously-submitted sample appears in the `seq_report` file, the flag `--update_reads` tells `ncbi_submit` to search for BioSasmple accessions of and include previously-submitted samples in the `submission.xml`. In most cases, if you are updating reads for a sample, a new SRA spuid is required. The `--spuid_endings` flag takes a parameter mapping samples that are being updated to a suffix. For any explicitely names samples, the suffix(es) will be added at the end of the automatically-generated SPUID. Usually '2' is a good suffix choice (unless another update has already been made using that same suffix for the sample of interest).

### Other metadata updates
These are not currently supported but could be added in the future if they seem important/useful.

##### Shell:
```bash
ncbi_submit file_prep \
    --config "${NCBI_CONFIG}" \
    --seq_report "${SEQ_REPORT}" \
    --outdir "${NCBI_DIR}" \
    --fastq_dir ${FASTQS} \
    --plate "${PLATE}" \
    --update_reads \
    --spuid_endings 'suffix1:samp1,samp2;suffix2:samp3'
```
##### Python:
```python
ncbi.write_presubmission_metadata(update_reads=True,spuid_endings={"sample1":"suffix1", "sample2":"suffix1", "sample3":"suffix2"})
```

***
## Input Paths
### Required Files:
  * `config`: Contains preset values and details about your lab, team, and submission plans that are necessary for submission.
  * `seq_report`: Main metadata file with sample details - can be equivalent to NCBI's BioSample TSV for use with the Submission Portal.
### Optional Files
  * `exclude_file`: Contains a list of "sample_name"s to exclude from NCBI submission (each one on a new line).
  * `barcode_map`: Used as a cross-reference. If all samples from `barcode_map` appear in `seq_report`, that's great. Otherwise, you'll get a warning with directions for adding samples to the `exclude_file` if they shouldn't be submitted. File should have no headers. Lines must be: "{barcode}\t{sample_name}".
### Sometimes Required Paths
  * `fastq_dir`: Required for `file_prep` and `ftp` if submitting reads to SRA. Indicates where the fastqs should be gathered from. Any fastqs with "sample_name" values that aren't supposed to be submitted will be ignored.
  * `outdir`: Highly recommended but will defualt to "./ncbi" or "./ncbi_test". A directory to house output (submission reports, `exclude_file`, output from `file_prep`). Will be created, if needed.
  * `subdir`: Only used for `ftp` tasks. A prefix to use for submissions for the given dataset. Defaults to `plate`, if plate is provided.

***
## Links to xml template examples/schema:
| File type | BioProject | BioSample | SRA | GenBank | Description/Link
|  --- | --- | --- | --- | --- | --- |
| Webpage | &check; | &check; | &check; | &check; | [Protocols & TSVs for use at Submission Portal](https://www.protocols.io/view/overview-of-ncbi-39-s-sars-cov-2-submission-proces-3byl476e2lo5/v5)
| XML | create | create | create |  | [SRA submission w/ new BioSample & BioProject](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/sra/samples/sra.submission.bs.bp.run.xml?view=co)
| XML | link | create | create |  | [SRA submission w/ new BioSample & existing BioProject](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/sra/samples/sra.submission.bs.run.xml?view=co)
| XML | link | link | create |  | [SRA submission w/ existing BioSample & BioProject](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/sra/samples/sra.submission.run.xml?view=co)
| XML |  |  |  | create | [GenBank XML](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/genbank/SARS-CoV-2/submission.xml?view=co)
| doc |  |  |  | example | [Example GenBank submission zip](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/genbank/SARS-CoV-2/)
| XSD |  | schema |  |  | [BioSample XML Schema](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/biosample/biosample.xsd?revision=71107&view=co)
| XSD | schema |  |  |  | [BioProject XML Schema](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/bioproject/bioproject.xsd?view=co)
| err | validate |  |  |  | [Submission Error Explanations](https://www.ncbi.nlm.nih.gov/projects/biosample/docs/submission/validation/errors.xml)
