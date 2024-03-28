#!/bin/bash

python ../../scripts/qpcr_specific_analysis.py -q input/qpcr_data_example.xlsx -qs Sheet1 -l input/qpcr_layout_example.xlsx -ls Sheet1 -f ../../artifacts/format_conversion.tsv -o output --pcopy 8.34e10 --fcopy 1.25e11