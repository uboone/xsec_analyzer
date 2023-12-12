FPM_Config="Configs/nuwro_file_properties.txt"
SYST_Config="Configs/systcalc.conf"
SLICE_Config="Configs/tutorial_slice_config.txt"
Univ_Output="/uboone/data/users/barrow/CC2P/Output.root"
PlotOutputDir="./Output/"

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
