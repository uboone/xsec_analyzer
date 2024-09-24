source setup_stv.sh

PROCESSED_NTUPLE_DIR="/exp/uboone/data/users/barrow/CC2P/"
UNIV_OUTPUT_FILE=${PROCESSED_NTUPLE_DIR}"/Universes.root"

MEASUREMENT_OUTPUT_FILE="./Output/"
UNF_MEAS_OUTPUT_FILE=${MEASUREMENT_OUTPUT_FILE}"/UnfoldedCrossSection.root"

PELEE_NTUPLE_CONFIG="./Configs/files_to_process.txt"
FPM_CONFIG="./Configs/file_properties.txt"
BIN_CONFIG="./Configs/tutorial_bin_config.txt"
XSEC_CONFIG="./Configs/xsec_config.txt"
SLICE_CONFIG="./Configs/tutorial_slice_config.txt"
SYST_CONFIG="./Configs/systcalc.conf"

./Scripts/ReprocessNTuples.sh ${PROCESSED_NTUPLE_DIR} ${PELEE_NTUPLE_CONFIG}
./Scripts/UniverseMaker.sh ${FPM_CONFIG} ${BIN_CONFIG} ${UNIV_OUTPUT_FILE} 
./Scripts/PlotSlices.sh ${FPM_CONFIG} ${SYST_CONFIG} ${SLICE_CONFIG} ${UNIV_OUTPUT_FILE} ${MEASUREMENT_OUTPUT_FILE}
./Scripts/Unfolder.sh ${XSEC_CONFIG} ${SLICE_CONFIG} ${UNF_MEAS_OUTPUT_FILE}
