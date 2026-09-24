# Raw-data utilities

Scripts in this folder are **not part of the main analysis pipeline** described
in the top-level [README.md](../README.md). They need the *raw* scanner data
(ismrmrd/Siemens raw MRI data with un-anonymized acquisition metadata), which
is not included in the OpenNeuro release ([ds004332](https://openneuro.org/datasets/ds004332))
that the rest of this repository works against.

## get_scan_end.py

Extracts scan start times from the raw ismrmrd headers, then computes scan
end times (start + a fixed per-sequence duration). **This is optional** --
the OpenNeuro release already ships the extracted scan end times
(`source/ScanEndTimes_*.txt`), so you only need this script if you want to
regenerate them from the raw data yourself.

### Getting the raw data

The raw data is hosted separately on PublicnEUro and is available **upon
request** (not open like the OpenNeuro release):

> PN000009 "Markerless Prospective Motion Correction"
> https://datacatalog.publicneuro.eu/dataset/PN000009%20Markerless%20Prospective%20Motion%20Correction/V1

### Requirements

- The Python [`ismrmrd`](https://pypi.org/project/ismrmrd/) package (not part of
  the main `requirements.txt`, since nothing else in this repo needs it):
  ```
  pip install ismrmrd
  ```
- Two environment variables:
  - `MOCO_DATASET_PATH` -- same as for the rest of the pipeline, pointing at
    the downloaded OpenNeuro ds004332 release. Used to look up each
    sequence's BIDS JSON sidecar (for e.g. checking which subjects have
    FLAIR/DWI/T2*) and as the base for the output directory
    (`$MOCO_DATASET_PATH/derivatives/results/Motion_Estimates/`).
  - `PUBLICNEURO_RAW_PATH` -- path to your local copy of the PublicnEUro raw
    data release (the directory that directly contains `sourcedata/`).

### Layout of the raw data

PublicnEUro is organized just like OpenNeuro: a `sourcedata/` folder holding
one subject directory per subject, each holding that subject's raw `.h5`
files directly (one file per acquisition), e.g.:

```
sourcedata/
  sub-001/
  sub-002/
  sub-003/
    meas_MID00132_FID18328_TCLmoco_off_still_T2_flair_sag_3D_tsevfl_harness_changed.h5
    meas_MID00019_FID15562_TCLmoco_off_nod_t1_mpr_3d_sag_p2_iso.h5
    meas_MID00059_FID15602_TCLmoco_on_still_ep2d_diff_exttracking.h5
    ...
```

The one difference from this repo's OpenNeuro side: subject IDs are
zero-padded to 3 digits (`sub-001/`) rather than 2 (`sub-01/`). Subject
numbering otherwise corresponds 1:1 between the two (confirmed) --
`get_scan_end.py` converts between the two paddings via the
`raw_to_bids_subject()` helper.

### Usage

```bash
export MOCO_DATASET_PATH=/path/to/ds004332-download/
export PUBLICNEURO_RAW_PATH=/path/to/publicneuro-download/
python3 get_scan_end.py
```

By default `get_scan_start = False` and `get_scan_end = True` inside the
script -- i.e. it expects `ScanTimes_*.txt` files (produced by a prior run
with `get_scan_start = True`) to already exist under
`$MOCO_DATASET_PATH/derivatives/results/Motion_Estimates/`, and only
computes the end times from them. Set `get_scan_start = True` for a full
from-scratch run against the raw data.
