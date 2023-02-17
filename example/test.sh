#!/usr/bin/env bash
main_dir=$(realpath $(dirname `dirname ${0}`))
plate="01-06-23-1-V2A"
plate_dir="${main_dir}/example/data/Clinical-01-06-23-1-V2A-fastqs"

# source login credentials
# note, the file can be created like this: `echo -e "ncbiUser='myUser'\nncbiPass='myPass'"" > .login_credentials`
. "${main_dir}/.login_credentials"
export ncbiUser=$ncbiUser
export ncbiPass=$ncbiPass

python "${main_dir}/ncbi_interact.py" \
    file_prep \
    --fastq_dir "${plate_dir}/fastqs/" \
    --seq_report "${plate_dir}/Sequencing_report-01-06-23-1-V2A-All.csv" \
    --plate "${plate}" \
    --outdir "${plate_dir}/ncbi/" \
    --primer_map "${plate_dir}/primer_map-01-06-23-1-V2A.csv" \
    --gisaid_log "${plate_dir}/gisaid_uploader.log" \
    # --primer_scheme 'vss2a'
