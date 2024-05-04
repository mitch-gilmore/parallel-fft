#!/bin/bash

# Define default remote path
DEFAULT_REMOTE_PATH="/projectnb/ec526/students/karabay"

# Check arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <local_path> [remote_path]"
    exit 1
fi

LOCAL_PATH=$1

# Check if remote path is provided, otherwise use default
if [ "$#" -eq 2 ]; then
    REMOTE_PATH=$2
else
    REMOTE_PATH=$DEFAULT_REMOTE_PATH
fi

REMOTE_USER="karabay"
REMOTE_HOST="scc1.bu.edu"

# Execute the scp command
scp -r $LOCAL_PATH $REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH

# Check if scp was successful
if [ $? -eq 0 ]; then
    echo "File/directory successfully transferred!"
else
    echo "Error occurred during transfer."
fi