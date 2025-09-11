#!/bin/bash
# Needs to be used inside the stage/bin folder 
PATH_EXEC=$1
OPENFACE_PATH=$2
VIDEO_PATH=$3
CSV_VIDEO_PATH=$4

set -m

# Added output redirection for stdout and stderr
${PATH_EXEC}action_unit_switch AUs/buste/ f001_head_color.tga ${OPENFACE_PATH} f001_head_normal.tga ${CSV_VIDEO_PATH} ${VIDEO_PATH} > filename 2>&1 &

# Wait for up to 5 seconds for the service to be ready.
for attempt in $(seq 1 20); do
    sleep 0.5
    if grep -q "Ready" filename; then
        echo "Launching Video"
        break
    fi
    if [[ attempt -eq 20 ]]; then
        echo "Error exec."
        exit
    fi
done
