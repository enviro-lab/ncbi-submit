#!/usr/bin/env bash
main_dir=$(realpath $(dirname `dirname ${0}`))
plate="01-06-23-1-V2A"
plate_dir="${main_dir}/example/data/Clinical-01-06-23-1-V2A-fastqs"

# source login credentials
# note, the file can be created with code like this: `echo -e "ncbiUser='myUser'\nncbiPass='myPass'"" > .login_credentials`
. "${main_dir}/.login_credentials"
export ncbiUser=$ncbiUser
export ncbiPass=$ncbiPass

### This script follows a typical submission order from start to finish
## Toggle `test_dir` and `test_mode` on/off from config, as needed

# Prepare files for XML submission
ncbi_submit \
    file_prep \
    --fastq_dir "${plate_dir}/fastqs/" \
    --seq_report "${plate_dir}/Sequencing_report-01-06-23-1-V2A-All.csv" \
    --plate "${plate}" \
    --outdir "${plate_dir}/ncbi/" \
    --primer_map "${plate_dir}/primer_map-01-06-23-1-V2A.csv" \
    --gisaid_log "${plate_dir}/gisaid_uploader.log"

# # Do XML submission for BioSample and SRA
# ncbi_submit \
#     ftp --submit --db 'bs_sra' \
#     --plate '01-06-23-1-V2A' \
#     --fastq_dir "${plate_dir}/fastqs/" \
#     --outdir "${plate_dir}/ncbi/" 

## Wait 10 minutes for first reports to become available
## Waiting longer may be necessary if BioSamples take longer to be created 
# sleep 600

# # Check submission and pull back reports (which contain submission status and BioSample accessions, if processed-ok)
# ncbi_submit \
#     ftp --check \
#     --plate '01-06-23-1-V2A' \
#     --fastq_dir "${plate_dir}/fastqs/" \
#     --outdir "${plate_dir}/ncbi/" 
#     # --db 'bs_sra' \

# # Add BioSample accessions to GenBank files
# ncbi_submit \
#     file_prep --prep_genbank \
#     --fasta "${plate_dir}/seqs-01-06-23-1-V2A.fasta" \
#     --fastq_dir "${plate_dir}/fastqs/" \
#     --seq_report "${plate_dir}/Sequencing_report-01-06-23-1-V2A-All.csv" \
#     --plate "${plate}" \
#     --outdir "${plate_dir}/ncbi/" \
#     --primer_map "${plate_dir}/primer_map-01-06-23-1-V2A.csv" \
#     --gisaid_log "${plate_dir}/gisaid_uploader.log" 
#     # --use_existing \
#     # --primer_scheme 'vss2a'

## Do XML submission for GenBank (Note: not tested)
# ncbi_submit \
#     ftp --submit --db 'gb' \
#     --plate '01-06-23-1-V2A' \
#     --fastq_dir "${plate_dir}/fastqs/" \
#     --outdir "${plate_dir}/ncbi/" 

# # Check on GenBank submission (Note: not tested)
# ncbi_submit \
#     ftp --check \
#     --db 'gb' \
#     --plate '01-06-23-1-V2A' \
#     --fastq_dir "${plate_dir}/fastqs/" \
#     --outdir "${plate_dir}/ncbi/" 