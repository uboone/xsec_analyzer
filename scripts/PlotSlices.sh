# Number of expected command-line arguments
num_expected=5

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: ./PlotSlices.sh FPM_Config SYST_Config SLICE_Config Univ_Output PlotOutputDir"
  exit 1
fi

FPM_Config=$1
SYST_Config=$2
SLICE_Config=$3
Univ_Output=$4
PlotOutputDir=$5

if [ ! -f "${FPM_Config}" ]; then
  echo "FPM_Config \"${FPM_Config}\" not found"
  exit 1
fi

if [ ! -f "${SYST_Config}" ]; then
  echo "SYST_Config \"${SYST_Config}\" not found"
  exit 2
fi

if [ ! -f "${SLICE_Config}" ]; then
  echo "SLICE_Config \"${SLICE_Config}\" not found"
  exit 3
fi

if [ ! -f "${Univ_Output}" ]; then
  echo "Universe File \"${Univ_Output}\" not found"
  exit 4
fi

if [ ! -d "${PlotOutputDir}" ]; then
  echo "Output directory \"${PlotOutputDir}\" not found"
  exit 5
fi

./Plotting/Slice_Plots ${FPM_Config} ${SYST_Config} ${SLICE_Config} ${Univ_Output} ${PlotOutputDir}
