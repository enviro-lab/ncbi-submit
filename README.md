# ncbi_interact
Submitting data to public databases is super important for publically funded laboratories, but it is not always a quick or intuitive process. NCBI Interact provides a simple and repeatable way to upload programmatic submissions to NCBI's SRA and GenBank with shared or unique BioProjects and BioSamples. Data can be uploaded as XML or zip files to either the Test or Production environments, and once there, the reports produced by NCBI can be analyzed to check on submission status and get BioSample accessions.

***
## Usage

`ncbi_interact.py` is intended for use on the command line, but the class `ncbi.NCBI` can be imported and used within custom python scripts.

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

### Setup
The required parameters vary by which of the above actions you're attempting but at minimum require a `plate` and `outdir`. To limit the number of parameters required via command line, a `config` file must be used When running from the command line, one of the three actions (`file_prep` or `ftp`) must be specified. With python, these are associated methods you may use on a single NCBI object.

#### Python instantiation (not needed on command line):
Note: This is the minimum required info for preparing data. Other parameters may be necessary for more functionality or other tasks.  
```python
from modules.ncbi import NCBI
ncbi = NCBI(
    fastq_dir = myFastqDir,
    seq_report = mySeqReport,
    plate = myPlate,
    outdir = myOutdir,
    config_file = myConfig,
    )
ncbi.write_presubmission_metadata()
```

### File Preperation
#### Shell:
```console
python ncbi_interact.py file_prep \
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
#### Python:
```python
ncbi.write_presubmission_metadata()
```

### File Submission
#### Shell:
```console
python ncbi_interact.py ftp \
    --submit --db bs_sra \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --outdir "${NCBI_DIR}" \
    -u "${ncbi_username}" \
    -p "${ncbi_password}" \
    --fastq_dir "${FASTQS}"
```
#### Python:
```python
ncbi.submit(db="bs_sra)
# wait awhile and try this to download reports and view submission status
ncbi.check(db="bs_sra)
```

### GenBank submission
(NOTE: not fully tested)
To link your fasta in GenBank to the associated reads, you'll want to add in the BioSample accessions before submitting. Follow these steps.
* Acquire BioSample accessions via one of these methods:
  * download accessions.tsv file from NCBI
    * (Do this if you submitted to BioSample via NCBI's Submission Portal)
  * use ncbi_interact.py
    * (Do this to avoid manual uploads  via NCBI's Submission Portal)
#### Shell:
```console
# dowload report.xml files to get accesssions from
python ncbi_interact.py ftp \
    --check --db ${DB} \
    --outdir "${NCBI_DIR}" \
    --config "${NCBI_CONFIG}" \
    -u "${ncbi_username}" \
    -p "${ncbi_password}" \
    --plate "${PLATE}" \
    --fastq_dir "${FASTQS}"

# add accessions to genbank.tsv
python ncbi_interact.py --prep_genbank \
    --outdir "${NCBI_DIR}" \
    --config ${NCBI_CONFIG} \
    --fasta "${GENERIC_CONSENSUS//PLATE/$PLATE}" \
    --plate "${PLATE}"

# submit to GenBank (NOTE: db='gb')
python ncbi_interact.py ftp \
    --submit --db gb \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --outdir "${NCBI_DIR}" \
    -u "${ncbi_username}" \
    -p "${ncbi_password}" \
    --fastq_dir "${FASTQS}"
```
#### Python:
```python
# dowload report.xml files to get accesssions from
ncbi.check(db="bs_sra)
# prepare genbank submission files and submit
ncbi.submit(db="gb")
# files can also be prepared without submitting via:
ncbi.write_genbank_submission_zip()
```

***
### Check Submission Status
Wait awhile (10+ minutes) for NCBI to start processing the submission. Then run this to download reports and view submission status.
This works for whichever db you want to check on. If not specified, you'll get results on all submitted dbs.

#### Shell:
```console
# check GenBank submission (NOTE: db='gb')
python ncbi_interact.py ftp \
    --check --db gb \
    --test_mode --test_dir \
    --config "${NCBI_CONFIG}" \
    --outdir "${NCBI_DIR}" \
    -u "${ncbi_username}" \
    -p "${ncbi_password}"
```
#### Python:
```python
# check GenBank submission (NOTE: db='gb')
ncbi.check(db="gb)
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
