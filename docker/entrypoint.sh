#!/bin/sh

eval "$(conda shell.bash hook)"
eval "conda activate mocohealthy"

eval "$@"
