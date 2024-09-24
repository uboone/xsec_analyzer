# Number of expected command-line arguments
num_expected=3

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: ./Unfolder.sh XSEC_CONFIG SLICE_CONFIG OUTPUT_FILE"
  exit 1
fi

XSEC_CONFIG=$1
SLICE_CONFIG=$2
OUTPUT_FILE=$3

if [ ! -f "$XSEC_CONFIG" ]; then
  echo "XSEC_CONFIG \"${XSEC_CONFIG}\" not found"
  exit 1
fi

if [ ! -f "$SLICE_CONFIG" ]; then
  echo "SLICE_CONFIG \"${SLICE_CONFIG}\" not found"
  exit 1
fi

Unfolder ${XSEC_CONFIG} ${SLICE_CONFIG} ${OUTPUT_FILE}
