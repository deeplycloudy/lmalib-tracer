# lmalib-tracer
Support code for using LMA data in the TRACER field campaign

## Installation

After cloning this repository, install with
``` conda env create environment.yml
conda activate lmatracer
pip install -e .
```

The realtime directory contains a script to sort and grid flashes that can be run with
```
python flash_sort.py LYLOUT_filenames*.dat.gz
```

If you use code in this repository in publications, please credit the authors:
- Eric Bruning (eric.bruning@ttu.edu)
