#!/bin/bash
# Run the Giotto-TDA pipeline and save the Nextflow log to the log/ folder.

mkdir -p log
nextflow -log "log/nextflow_$(date +%Y-%m-%d_%H-%M-%S).log" run main.nf "$@"
