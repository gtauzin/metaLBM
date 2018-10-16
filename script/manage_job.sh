#!/bin/bash

job_type=$1
hostname=$2
input_json=$3

job_name_list=job_name_${hostname}.csv

job_id=$(cat "${job_name_list}" | wc -l)

job_name=${job_type}-${job_id}

echo "${input_json}, ${job_name}" >> ${job_name_list}

echo $job_name
