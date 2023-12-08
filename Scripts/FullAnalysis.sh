source setup_stv.sh

PROCESSED_NTUPLE_DIR="/uboone/data/users/barrow/CC2P/"
UNIV_OUTPUT_FILE=${PROCESSED_NTUPLE_DIR}"/Output.root"

PELEE_NTUPLE_CONFIG="./Configs/files_to_process.txt"
FPM_CONFIG="./Configs/file_properties.txt"
BIN_CONFIG="./Configs/tutorial_bin_config.txt"

./Scripts/ReprocessNTuples.sh ${PROCESSED_NTUPLE_DIR} ${PELEE_NTUPLE_CONFIG}
./Scripts/UniverseMaker.sh ${FPM_CONFIG} ${BIN_CONFIG} ${UNIV_OUTPUT_FILE}

