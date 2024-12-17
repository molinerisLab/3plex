#!/bin/bash


# Default values
TARGET_PATH="3plex.sif"
MOUNT_DIRECTORY=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --TARGET_PATH) TARGET_PATH="$2"; shift ;;
        --MOUNT_DIRECTORY) MOUNT_DIRECTORY="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done


# If MOUNT_DIRECTORY is not empty, mount it
if [ -n "$MOUNT_DIRECTORY" ]; then
    singularity shell -B ./:/3plex -B $MOUNT_DIRECTORY:$MOUNT_DIRECTORY --pwd /3plex "$TARGET_PATH"/ 
else
   singularity shell -B ./:/3plex --pwd /3plex "$TARGET_PATH"/ 
fi
