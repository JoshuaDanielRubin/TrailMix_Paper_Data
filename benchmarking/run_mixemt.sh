#!/bin/bash

# Get the input file name
input_file="$1"
output_prefix="$2"
# Define the log file
log_file="${output_prefix}.log"
# Run the mixemt command with the required options and redirect both stdout and stderr to the log file
python mixemt/bin/mixemt "$input_file" --ref rCRS.fa --seed 42 -v -V -s "$output_prefix"
