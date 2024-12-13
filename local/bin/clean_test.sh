#!/bin/bash

# Check if the target name is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <target>"
  exit 1
fi

TARGET=$1

# Remove all output files listed in the summary
snakemake $TARGET --summary | tail -n +2 | awk '{print $1}' | xargs rm -f

# Remove log files
find . -name "*.log" -type f -delete

# Remove empty directories
find . -type d -empty -delete

# Remove pato_e.txt
rm pato_e.txt