#!/bin/bash

# Check if a command-line argument is provided
if [ -z "$1" ]; then
    # No argument provided, use default path
    TARGET_PATH="3plex.sif"
else
    # Use the provided argument as the target path
    TARGET_PATH="$1"
fi

# Run the singularity build command
singularity build "$TARGET_PATH" .singularity/singularity.def