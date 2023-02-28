#!/usr/bin/env python
# NOTE: This has been adapted from the config used by MultiSub
# Multisub allows for submission to GISAID, NCBI, and something else
# NCBI-Interact goes more in depth into NCBI programmatic submissions.

# For testing purposes (both default to false if commented out) 
test_dir = True # Prepare data with test accession and submit to NCBI's Test directory
# test_mode = True # Do everything except the actual upload of data

# --- For automated NCBI uploads via FTP --
# Some of the below values will be (not yet) used to populate the NCBI submission template file
# create your own at: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
# `host` defaults to "ftp-private.ncbi.nlm.nih.gov", but can be adjusted if submitting elsewhere
template = "./template.sbt"
host = ""
ncbiUser = ""
ncbiPass = ""
centerAbbr = "CLT"

# extra variables that can be arguments instead
controls = "NTC|Neg|Pos"

# Filename used by ncbi_prep.py to check for samples that have already been submitted. 
# Any sample names found that have already been submitted will cause the script to fail.
# If that sample should be ignored, it can be added to the file "`outdir`/samples2exclude.txt"
# Or,  if you're sure it needs to be submitted afresh, remove it from `submitted_samples`
# Default: unless a filename is provided, the script will not check for previous submissions.
# submitted_samples = ""

# general info
# format: (lastname, firstname, middle initial)
contact = ("Researcher", "Agood", "")
email = "aresearch@uni.edu"
phone = "123-456-7890"

# full name of country, used for NCBI export, e.g. France
# see https://www.ncbi.nlm.nih.gov/genbank/collab/country/
ncbiCountry = "USA"
# alternative email address, only used for NCBI submissions and not required
# altEmail = "johndoe@gmail.com"

affiliation = {
"affil" : "University of Education",
"div" : "University of Education Environmental Monitoring Laboratory",
"city" : "Somewhere",
"sub" : "NC",
"country" : "USA", 
"street" : "123 Any Street",
"email" : email,
"phone" : phone,
"postal-code" : "24680",
}

authors = [
    # format: (lastname, firstname, middle initial)
    ("Researcher", "Agood", ""),
    ("Seymore", "Example", ""),
]

# --- REQUIRED CONFIGURATION ENDS ---
# below you find configuration sections for each export database

# -- Defaults for NCBI submissions ---
# NOTE: genbank metadata will be compiled from values in biosample/sra attributes
# NOTE: empty string values override default derived values from the seq_report csv
# NOTE: any commented out values will not be included in the TSV or XML files
# NOTE: any empty strings will not be included in the XML files but will appear empty or as 'nan' in the TSVs

bioproject = {}
# A BioProject is an entry for your project, it has a title, authors, possibly papers and groups together all your submissions
# Instructions on how to create one are here: https://www.protocols.io/view/sars-cov-2-ncbi-submission-protocol-sra-biosample-bui7nuhn?step=3
# A new one will be created if no BioProject accession exists in the config
#   or you can create one here https://submit.ncbi.nlm.nih.gov/subs/bioproject/
# REQUIRED BioProject attributes
bioproject['bioproject_accession'] = "PRJNA111111"
bioproject['bioproject_test_accession'] = "PRJNA553747"
bioproject['create_new'] = False
# OPTIONAL BioProject attributes - use if creating a new BioSample
bioproject['spuid'] = "SARS-CoV-2 WGS: University of Education Environmental Monitoring Laboratory"
    # spuid must be unique to submission center
bioproject["title"] = "Genome sequencing of SARS-CoV-2 as part of the University of Education campus COVID-19 wastewater surveillance effort"
bioproject["description"] = "Genome sequencing of SARS-CoV-2 as part of the University of Education campus COVID-19 surveillance effort"
    # Very short descriptive name of the project for caption, labels, etc. For example: 1000 Genomes Project
bioproject["website_name"] = ""
bioproject["url"] = ""
bioproject["dataType"] = "genome sequencing and assembly" # if assembling
# bioproject["dataType"] = "genome sequencing" # if not assembling
bioproject["organism"] = "Severe acute respiratory syndrome coronavirus 2"
bioproject["scope"] = "eMultiisolate"


biosample = {}
# REQUIRED BioSample attributes:
# Pick the correct package and organism for your data type: https://www.ncbi.nlm.nih.gov/biosample/docs/packages/. Here are the SARS-CoV-2 options
biosample['package'] = "SARS-CoV-2.cl.1.0"
biosample['organism'] = "Severe acute respiratory syndrome coronavirus 2"
# biosample['package'] = "SARS-CoV-2.wwsurv.1.0"
# biosample["organism"] = "wastewater metagenome"
# The mandatory attributes are listed here: https://www.ncbi.nlm.nih.gov/biosample/docs/packages/SARS-CoV-2.cl.1.0/
biosample['collected_by'] = affiliation["div"]
biosample['geo_loc_name'] = "USA: North Carolina"
biosample['host'] = "Homo sapiens"
biosample['host_disease'] = "COVID-19"
biosample['isolation_source'] = "Clinical"
biosample['bioproject_accession'] = bioproject['bioproject_accession']

# OPTIONAL BioSample attributes:
# list (comma-delimited) of desired fields to go in BioSample submission. This is useful if --seq_report is already formatted for BioSample submission: just specify all fields present there. (empty fields will be dropped)
# biosample['all_cols'] = 'field1,field2,..'
# NOTE: more of these can be added just by changing the number at the end. Using these plus biosample['all_cols'] might cause issues. Best to use one or the other
biosample['sars_cov_2_diag_gene_name_1'] = 'N (orf9)' #"N Gene"
# biosample['sars_cov_2_diag_gene_name_2'] = ''
# Corresponding number and spreasheet column name for the ct value
# as dict: {1:"val1",2:"val2",...}
sars_cov_2_diag_pcr_ct_values = {1:'Ct N gene'}

# spreadsheet_col__sars_cov_2_diag_pcr_ct_value_2 = ''
biosample['sample_title'] = '' # empty string values override default derived values from the seq_report csv # NOTE: we don't want sample_title to show because that ID is private for us
biosample['antiviral_treatment_agent'] = "Restricted Access"
biosample['collection_device'] = "Swab"
biosample['collection_method'] = "Swabbing"
biosample['date_of_prior_antiviral_treat'] = "Restricted Access"
biosample['date_of_prior_sars_cov_2_infection'] = "Restricted Access"
biosample['date_of_sars_cov_2_vaccination'] = "Restricted Access"
biosample['exposure_event'] = "Restricted Access"
biosample['geo_loc_exposure'] = "Restricted Access"
biosample['host_age'] = "Restricted Access"
biosample['host_anatomical_material'] = "Tissue"
biosample['host_anatomical_part'] = "Nasopharynx"
biosample['host_body_product'] = "Not Applicable"
biosample['host_disease_outcome'] = "Restricted Access"
biosample['host_health_state'] = "Restricted Access"
biosample['host_recent_travel_loc'] = "Not Collected"
biosample['host_recent_travel_return_date'] = "Not Collected"
biosample['host_sex'] = "Missing"
biosample['host_specimen_voucher'] = "Not Collected"
biosample['host_subject_id'] = "Restricted Access"
biosample['lat_lon'] = "Missing"
biosample['passage_method'] = "Not Collected"
biosample['passage_number'] = "Not Collected"
biosample['prior_sars_cov_2_antiviral_treat'] = ""
biosample['prior_sars_cov_2_infection'] = ""
biosample['prior_sars_cov_2_vaccination'] = ""
biosample['purpose_of_sampling'] = "Diagnostic Testing"
biosample['purpose_of_sequencing'] = "Targeted Surveillance"
biosample['sequenced_by'] = "University of Education Environmental Monitoring Laboratory"
biosample['vaccine_received'] = "Restricted Access"
biosample['virus_isolate_of_prior_infection'] = "Restricted Access"
biosample['description'] = ""

genbank = {}
# REQUIRED GenBank presets:
genbank["comment"] = {
    # indicates whether genbank CSV & XML should be skipped (True) or produced (False)
    "include_comment":True,
    "data":{}
}
# If including the comment file add any fields and preset values here. Not set up for individualized differences between samples
genbank["comment"]["data"]["sequencing_technology"] = "PromethION"
genbank["comment"]["data"]["assembly_method"] = "ARTIC-nCoV-bioinformaticsSOP"
genbank["comment"]["data"]["assembly_version_or_date"] = "1.1.0"

sra = {}
# REQUIRED SRA attributes:
sra['title'] = "Genomic sequencing of SARS-CoV-2: Clinical, USA:North Carolina"
sra['library_strategy'] = 'AMPLICON'
sra['library_source'] = 'VIRAL RNA'
sra['library_selection'] = 'PCR'
sra['library_layout'] = 'single'
sra['platform'] = 'OXFORD_NANOPORE'
sra['instrument_model'] = 'PromethION'
sra['filetype'] = 'fastq'
# OPTIONAL SRA attributes
sra['sequence_submitter_contact_email'] = email
sra['amplicon_size'] = ''
sra['raw_sequence_data_processing_method'] = ''
sra['dehosting_method'] = 'kraken2'
sra['loader'] = 'fastq-load.py' # this is used specifically for ONT sequencing

## ----- Variable Primer Scheme and Protocol Info ------
allowed_schemes = ["V3","V3.1","V4","V4.1","vsl1a","vss1a","vss2a","vss2b"]
# all below schemes must be in the allowed_schemes list to prevent errors
protocol_scheme={scheme:{} for scheme in allowed_schemes}

scheme="V3"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="https://github.com/joshquick/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/"
protocol_scheme[scheme]["design_description"]= "Viral sequencing was performed following a tiling amplicon strategy using the V3 primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies. Libraries were prepared using High-Capacity cDNA Reverse Transcription Kit by ThermoFisher, Q5® High-Fidelity DNA Polymerase pack by NEB, EB Buffer by Omega Biotek, NEBNext® Ultra™ II End Repair/dA-Tailing Module ,NEBNext® Ultra™ II Ligation Module and NEBNext® Quick Ligation Module by NEB, Sequencing auxiliary vials by ONT, Adapter Mix II Expansion by ONT."
protocol_scheme[scheme]["sequencing_protocol_name"]= "A SARS-CoV-2 Surveillance Sequencing Protocol Optimized for Oxford Nanopore PromethION"

scheme="V3.1"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="Combination of https://github.com/joshquick/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/ and the patch added to V4 by https://github.com/joshquick/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4.1"
protocol_scheme[scheme]["design_description"]= "Viral sequencing was performed following a tiling amplicon strategy using the V3.1 primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies. Libraries were prepared using High-Capacity cDNA Reverse Transcription Kit by ThermoFisher, Q5® High-Fidelity DNA Polymerase pack by NEB, EB Buffer by Omega Biotek, NEBNext® Ultra™ II End Repair/dA-Tailing Module ,NEBNext® Ultra™ II Ligation Module and NEBNext® Quick Ligation Module by NEB, Sequencing auxiliary vials by ONT and Adapter Mix II Expansion by ONT."
protocol_scheme[scheme]["sequencing_protocol_name"]="A SARS-CoV-2 Surveillance Sequencing Protocol Optimized for Oxford Nanopore PromethION"

scheme="V4"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="https://github.com/joshquick/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4"
protocol_scheme[scheme]["design_description"]= "Viral sequencing was performed following a tiling amplicon strategy using the V4 primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies. Libraries were prepared using High-Capacity cDNA Reverse Transcription Kit by ThermoFisher, Q5® High-Fidelity DNA Polymerase pack by NEB, EB Buffer by Omega Biotek, NEBNext® Ultra™ II End Repair/dA-Tailing Module ,NEBNext® Ultra™ II Ligation Module and NEBNext® Quick Ligation Module by NEB, Sequencing auxiliary vials by ONT, Adapter Mix II Expansion by ONT."
protocol_scheme[scheme]["sequencing_protocol_name"]= "A SARS-CoV-2 Surveillance Sequencing Protocol Optimized for Oxford Nanopore PromethION"

scheme="V4.1"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="https://github.com/joshquick/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4.1"
protocol_scheme[scheme]["design_description"]= "Viral sequencing was performed following a tiling amplicon strategy using the V4.1 primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies. Libraries were prepared using High-Capacity cDNA Reverse Transcription Kit by ThermoFisher, Q5® High-Fidelity DNA Polymerase pack by NEB, EB Buffer by Omega Biotek, NEBNext® Ultra™ II End Repair/dA-Tailing Module ,NEBNext® Ultra™ II Ligation Module and NEBNext® Quick Ligation Module by NEB, Sequencing auxiliary vials by ONT, Adapter Mix II Expansion by ONT."
protocol_scheme[scheme]["sequencing_protocol_name"]= "A SARS-CoV-2 Surveillance Sequencing Protocol Optimized for Oxford Nanopore PromethION"

scheme="vsl1a"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="https://github.com/nebiolabs/VarSkip/tree/main/schemes/NEB_VarSkip/V1a-long"
protocol_scheme[scheme]["design_description"]="Viral sequencing was performed following a tiling amplicon strategy using the vsl1a primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies . Libraries were prepared using NEBNext® ARTIC SARS-CoV-2 Companion Kit"
protocol_scheme[scheme]["sequencing_protocol_name"]= "Instruction manual:NEBNext® ARTIC SARS-CoV-2 Companion Kit: Chapter 3 along with Clean up after PCR amplification from Chapter 2"

scheme="vss1a"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="https://github.com/nebiolabs/VarSkip/tree/main/schemes/NEB_VarSkip/V1a"
protocol_scheme[scheme]["design_description"]="Viral sequencing was performed following a tiling amplicon strategy using the vss1a primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies . Libraries were prepared using NEBNext® ARTIC SARS-CoV-2 Companion Kit"
protocol_scheme[scheme]["sequencing_protocol_name"]= "Instruction manual:NEBNext® ARTIC SARS-CoV-2 Companion Kit: Chapter 3 along with Clean up after PCR amplification from Chapter 2"

scheme="vss2a"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="https://github.com/nebiolabs/VarSkip/tree/main/schemes/NEB_VarSkip/V2a"
protocol_scheme[scheme]["design_description"]="Viral sequencing was performed following a tiling amplicon strategy using the vss2a primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies . Libraries were prepared using NEBNext® ARTIC SARS-CoV-2 Companion Kit"
protocol_scheme[scheme]["sequencing_protocol_name"]= "Instruction manual:NEBNext® ARTIC SARS-CoV-2 Companion Kit: Chapter 3 along with Clean up after PCR amplification from Chapter 2"

scheme="vss2b"
protocol_scheme[scheme]["amplicon_PCR_primer_scheme"]="https://github.com/nebiolabs/VarSkip/tree/main/schemes/NEB_VarSkip/V2b"
protocol_scheme[scheme]["design_description"]="Viral sequencing was performed following a tiling amplicon strategy using the vss2b primer scheme. Sequencing was performed using a PromethIon by Oxford Nanopore Technologies . Libraries were prepared using NEBNext® ARTIC SARS-CoV-2 Companion Kit"
protocol_scheme[scheme]["sequencing_protocol_name"]= "Instruction manual:NEBNext® ARTIC SARS-CoV-2 Companion Kit: Chapter 3 along with Clean up after PCR amplification from Chapter 2"

# if uncommented, these will override any values specified in the mapping above
# sra["amplicon_PCR_primer_scheme"] = ""
# sra["design_description"] = ""
# sra["sequencing_protocol_name"] = ""
