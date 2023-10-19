#!/bin/bash

# Define the full path to the uvicorn executable
uvicorn_executable="/home/ubuntu/.pyenv/shims/uvicorn"
cd /home/ubuntu/plasmasphere/DC-BPIM/API
"$uvicorn_executable" apis:app --reload --port 8080 --host 0.0.0.0
