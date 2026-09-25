#!/bin/bash
#
# Runs the complete analysis pipeline end to end against a downloaded
# OpenNeuro ds004332 release: sets up the Python environment, then image
# quality assessment, motion analysis, cortical thickness (Fig. 8 -- by
# far the slowest step, hours to days depending on cores), and all
# remaining manuscript plots (Figs. 2, 4-7 and the Figure-3-style example).
#
# Usage:
#   bash run_full_analysis.sh <ds004332-download-path> <freesurfer-home> [analysis_cort_thickness.py args...]
#
# or, with the environment variables already exported:
#   export MOCO_DATASET_PATH=/path/to/ds004332-download/
#   export FREESURFER_HOME=/path/to/freesurfer
#   bash run_full_analysis.sh [analysis_cort_thickness.py args...]
#
# Any trailing arguments are passed through to analysis_cort_thickness.py,
# e.g. to control parallelism or run on a subject subset:
#   bash run_full_analysis.sh /data/ds004332-download /opt/freesurfer --jobs 8 --threads-per-job 2
#   bash run_full_analysis.sh /data/ds004332-download /opt/freesurfer --subjects sub-01 sub-02 sub-03
#
# Note: do NOT `set -u` in this script or anything it sources -- FreeSurfer's
# own SetUpFreeSurfer.sh references an unbound variable and breaks under it.
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Positional args (if they're existing directories) set the env vars.
if [ $# -ge 1 ] && [ -d "$1" ]; then
    export MOCO_DATASET_PATH="$1"
    shift
fi
if [ $# -ge 1 ] && [ -d "$1" ]; then
    export FREESURFER_HOME="$1"
    shift
fi
# Anything left in "$@" is passed through to analysis_cort_thickness.py.

if [ -z "${MOCO_DATASET_PATH:-}" ]; then
    echo "Usage: $0 <ds004332-download-path> <freesurfer-home> [analysis_cort_thickness.py args...]"
    echo "  (or export MOCO_DATASET_PATH and FREESURFER_HOME yourself first)"
    exit 1
fi
if [ -z "${FREESURFER_HOME:-}" ]; then
    echo "FREESURFER_HOME is not set (pass it as the 2nd argument, or export it)."
    exit 1
fi
case "$MOCO_DATASET_PATH" in
    */) ;;
    *) export MOCO_DATASET_PATH="${MOCO_DATASET_PATH}/" ;;
esac

echo "MOCO_DATASET_PATH=$MOCO_DATASET_PATH"
echo "FREESURFER_HOME=$FREESURFER_HOME"
echo "analysis_cort_thickness.py extra args: $*"

export MPLBACKEND=Agg
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"

if [ ! -d .venv ]; then
    echo "Creating Python venv at $SCRIPT_DIR/.venv ..."
    python3 -m venv .venv
fi
source .venv/bin/activate
pip install -q -r requirements.txt

run_step () {
    echo ""
    echo "=================================================================="
    echo "STEP: $* ($(date))"
    echo "=================================================================="
    "$@"
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "STEP FAILED (exit $rc): $*"
        exit $rc
    fi
}

run_step python3 analysis_img_quality.py
run_step python3 analysis_motion_data.py
run_step python3 analysis_cort_thickness.py "$@"
run_step python3 plot_generation_rewrite.py
run_step python3 make_figure3_example.py

echo ""
echo "=================================================================="
echo "Full analysis complete ($(date))."
echo "Results under ${MOCO_DATASET_PATH}derivatives/results/"
echo "=================================================================="
