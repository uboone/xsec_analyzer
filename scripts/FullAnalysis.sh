source setup_stv.sh

PROCESSED_NTUPLE_DIR="/exp/uboone/data/users/barrow/CC2P/"
UNIV_OUTPUT_FILE=${PROCESSED_NTUPLE_DIR}"/Universes.root"

MEASUREMENT_OUTPUT_FILE="./Output/"
UNF_MEAS_OUTPUT_FILE=${MEASUREMENT_OUTPUT_FILE}"/UnfoldedCrossSection.root"

PELEE_NTUPLE_CONFIG="./configs/files_to_process.txt"
FPM_CONFIG="./configs/file_properties.txt"
BIN_CONFIG="./configs/tutorial_bin_config.txt"
XSEC_CONFIG="./configs/xsec_config.txt"
SLICE_CONFIG="./configs/tutorial_slice_config.txt"
SYST_CONFIG="./configs/systcalc.conf"

./Scripts/ReprocessNTuples.sh ${PROCESSED_NTUPLE_DIR} ${PELEE_NTUPLE_CONFIG}
./Scripts/UniverseMaker.sh ${FPM_CONFIG} ${BIN_CONFIG} ${UNIV_OUTPUT_FILE} 
./Scripts/PlotSlices.sh ${FPM_CONFIG} ${SYST_CONFIG} ${SLICE_CONFIG} ${UNIV_OUTPUT_FILE} ${MEASUREMENT_OUTPUT_FILE}
./Scripts/Unfolder.sh ${XSEC_CONFIG} ${SLICE_CONFIG} ${UNF_MEAS_OUTPUT_FILE}
