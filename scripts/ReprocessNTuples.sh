#!/bin/bash

# Number of expected command-line arguments
num_expected=2

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: ./ReprocessNTuples.sh OUTPUT_DIRECTORY NTUPLE_LIST_FILE"
  exit 1
fi

output_dir=$1
ntuple_list_file=$2

if [ ! -f "$ntuple_list_file" ]; then
  echo "Ntuple list file \"${ntuple_list_file}\" not found"
  exit 1
fi

if [ ! -d "${output_dir}" ]; then
  echo "Output directory \"${output_dir}\" not found"
  exit 2
fi

input_files=()
# Loop over each line of the ntuple list file
while read line; do
  # Select lines that do not begin with a '#' character and contain at least
  # one non-whitespace character. These are assumed to be input file names
  if [[ ! $line = \#* ]] && [[ $line = *[^[:space:]]* ]]; then
    # Process the next input ntuple file 
      input_files+=(${line})
  fi
done < "${ntuple_list_file}"

# Calculate total number of input files
total_files=${#input_files[*]}
echo "Total number of files = "${total_files}

counter=0
# Loop over each input file
for file in "${input_files[@]}"
do
    input_file_name=${file}
    output_file_name="${output_dir}/xsec-ana-$(basename ${input_file_name})"
    echo "Starting file:"${counter}"/"${total_files}
    date
    time ProcessNTuples ${input_file_name} ${output_file_name}
    date
    counter=$((counter + 1))
done
