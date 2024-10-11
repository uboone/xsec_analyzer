#!/bin/bash

# Function to print usage information
usage() {
  echo "Usage: $0 [-o <out directory>] [-v version] [-r runnumbers] [-s samples]"
  echo "  -o directory   Specify the output directory."
  echo "  -v version     version of your reprocess ntuples"
  echo "  -r runnumbers  Specify runnumbers to filter (comma-separated, e.g., 1,2,3)."
  echo "  -s samples     Specify runnumbers to filter "
  echo "                      run1 : numuMC,nueMC,dirtMC,extBNB,onBNB,openBNB,detVarCV,detVarLYDown,detVarLYRayleigh,detVarLYAttenuation,detVarSEC,detVarRecomb2,detVarModX,detVarModYZ,"
  echo "  -v version     Specify the version in teck note (e.g. v00_00_01)"
  exit 1
}

version=current
config_dir=${XSEC_ANALYZER_DIR}/configs
input_file=${config_dir}/files_to_process_cthorpe.txt
output_config=file_properties_cthorpe
output_dir=../Output
filter_runnumbers=1,2,3,4,5
filter_samples="numuMC,nueMC,dirtMC,extBNB,onBNB,openBNB,detVarCV,detVarLYdown,detVarLYrayl,detVarLYatten,detVarSCE,detVarRecomb2,detVarWMX,detVarWMYZ,detVarWMAngleXZ,detVarWMAngleYZ,detVarCVExtra,altCVMC"
# Parse command-line arguments
while getopts ":o:r:s:v:" opt; do
  case ${opt} in
    o )
      output_dir=$OPTARG
      ;;
    r )
      filter_runnumbers=$OPTARG
      ;;
    s )
      filter_samples=$OPTARG
      ;;
    v )
      version=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

# Function to check the existence of a file or directory
check_existence() {
    if [ -z "$1" ]; then
        echo "Usage: check_existence <path>"
        return 1
    fi

    if [ -e "$1" ]; then
        if [ -f "$1" ]; then
            echo "The file '$1' exists. Please remove it or update the version."
            exit 1
        elif [ -d "$1" ]; then
            echo "The path '$1' exists. Please remove it or select a new path."
            exit 1
        fi
    else
        echo "The path '$1' does not exist. Continue..."
    fi
}

check_existence_dir(){

  if [ -z "$1" ]; then
    echo "Usage: check_existence <path>"
    return 1
  fi

    if [ -e "$1" ]; then
        if [ -d "$1" ]; then
            echo "The path '$1' exists. Please remove it or select a new path."
            exit 1
        fi
    else
        echo "The path '$1' does not exist. Create it..."
        mkdir -p $1
    fi
}


# Check if input file is provided
if [ -z "$input_file" ]; then
  usage
fi

# Convert comma-separated runnumbers to array
if [ -n "$filter_runnumbers" ]; then
  IFS=',' read -r -a runnumbers <<< "$filter_runnumbers"
fi
# Convert comma-separated runnumbers to array
if [ -n "$filter_runnumbers" ]; then
  IFS=',' read -r -a samples <<< "$filter_samples"
fi

# Function to check if a value is in an array
contains() {
  local array="$1[@]"
  local seeking=$2
  local in=1
  for element in "${!array}"; do
    if [[ $element == "$seeking" ]]; then
      in=0
      break
    fi
  done
  return $in
}

# Read and parse the file



output_dir=`realpath ${output_dir}`
output_config=${output_config}_${version}.txt

#check_existence ${config_dir}/${output_config}
#check_existence_dir ${output_dir}


if [ ! -f "$input_file" ]; then
  echo "Ntuple list file \"${input_file}\" not found"
  exit 1
fi

if [ ! -d "${output_dir}" ]; then
  echo "Output directory \"${output_dir}\" not found"
  exit 2
fi


echo "# `date`" > ${config_dir}/${output_config}
while IFS= read -r line; do
  # Skip comments and empty lines
  if [[ "$line" =~ ^# ]] || [[ -z "$line" ]]; then
  echo $line >> ${config_dir}/${output_config}
    continue
  fi
  
  # Split the line into parts
  parts=($line)
  file_name=${parts[0]}
  index=${parts[1]}
  sample_type=${parts[2]}
  num_events=${parts[3]:-}
  scaling_factor=${parts[4]:-}
  
  # Filter by runnumbers if specified
  if [ -n "$filter_runnumbers" ]; then
    if ! contains runnumbers "$index"; then
      continue
    fi
  fi
 
  # Filter by runnumbers if specified
# echo $sample_type
  if [ -n "$filter_runnumbers" ]; then
    if ! contains samples "$sample_type"; then
      continue
    fi
  fi

  echo  -n "$sample_type,"
  # Print the parsed information
  echo "File Name: $file_name"
  echo "Index: $index"
  echo "Sample Type: $sample_type"
  echo "Number of Events: $num_events"
  echo "Scaling Factor: $scaling_factor"
  output_file_name="${output_dir}/xsec-${version}-$(basename ${file_name})"
  echo "Output file name: ${output_file_name}"
  echo "time ProcessNTuples ${file_name} JOINTCC0pi ${output_file_name}"
  time ProcessNTuples ${file_name} JOINTCC0pi ${output_file_name}
  echo "${output_file_name} $index $sample_type $num_events $scaling_factor" >> ${config_dir}/${output_config}
  echo "--------------------------------"

done < "$input_file"
echo 
