#!/bin/bash

delta=$1

sed -i "s/delta\s*=\s*[0-9]*\.[0-9]*/delta = $delta/" FCI.py
sed -i "s/vib_run_dir=delta_[0-9]*\.[0-9]*/vib_run_dir=delta_$delta/" directory_setup.sh
sed -i "s/delta=[0-9]*\.[0-9]*/delta=$delta/" run_vib.sh




