python -m venv env
source env/bin/activate
which python
pip install --upgrade pip
pip install wheel setuptools
pip install -r requirements.txt

export PYTHONPATH=$PYTHONPATH:$PWD
export MPLCONFIGDIR=$PWD/.config/matplotlib