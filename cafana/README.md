## CAFAna

CAFAna is a framework for parsing the standard output of SBN CAF files. [sbnana](https://github.com/SBNSoftware/sbnana) is the main repository for the SBN software suite, which contains `cafe`. 

## Setup
Setup `sbnana` to get access to `cafe` executable.
```bash
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
ups list -aK+ sbnana #Use most recent version > 10_00_00 
setup sbnana v10_00_00 -q e26:prof
```

## Usage
The macros should be setup to store the event classifications into `event_type` branch. `is_signal` identifies if the event is reconstructed as signal. These are both essential for `ProcessNTuples` to work.
```bash
cafe -bq <macro.C>
```

## Example

```bash
cafe -bq xsec_analyzer.C
```

## Output

The output is a ROOT file with the name `xsec_analyzer.root`. The file contains a tree with the name `xsec_analyzer`. This is used as input for `ProcessNTuples` to feed into downstream fitting tools.
