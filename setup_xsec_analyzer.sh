#!/bin/bash

# Get the name of the OS release (used to distinguish between SL7 and EL9)
if test -f "/etc/os-release"; then
  MY_OS_REL=$(cat /etc/os-release | grep ^NAME \
    | sed -e 's/NAME=//g' -e 's/"//g')
fi

# Sets up the local environment for working with xsec_analyzer
#if [ "$MY_OS_REL" = "AlmaLinux" ]; then
  # On AL9, we set up ROOT and a recent compiler version
#  source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
#  spack load gcc@12.2.0 arch=linux-almalinux9-x86_64_v3
#  spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
#elif [ "$MY_OS_REL" = "Scientific Linux" ]; then
  # On SL7, we get ROOT as a side-effect of setting up uboonecode
  source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
  setup uboonecode v08_00_00_84 -q e17:prof
elif [ "$(uname)" == "Darwin" ]; then
  # We are on macOS, so tell the user they're on their own
  echo "macOS system detected, unable to automatically set up ROOT"
else
  echo "WARNING: Unrecognized OS name \"${MY_OS_REL}\""
  echo "Unable to automatically set up ROOT"
fi

# Finds the directory where this script is located. This method isn't
# foolproof. See https://stackoverflow.com/a/246128/4081973 if you need
# something more robust for edge cases (e.g., you're calling the script using
# symlinks).
THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export XSEC_ANALYZER_DIR=${THIS_DIRECTORY}
export PATH=${PATH}:${XSEC_ANALYZER_DIR}/bin:${HOME}/.local/bin
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${XSEC_ANALYZER_DIR}/include

echo "$(uname)"

# Set the library path for loading the XSecAnalyzer shared library at runtime
if [ "$(uname)" = "Darwin" ]; then
  # macOS platform (use both variables just in case)
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THIS_DIRECTORY}/lib
  export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THIS_DIRECTORY}/lib
else
  # Assume a GNU/Linux platform
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THIS_DIRECTORY}/lib
fi
