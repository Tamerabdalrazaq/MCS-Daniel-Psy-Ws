# Do previous target locations bias our attention during search in natural scenes? reanalysis of the MCS dataset

This repository contains the analysis code and processed data for our two eye-tracking visual search studies:

- **Study 1 – Microwaves** (`m`)
- **Study 2 – Clocks** (`c`)

---

## 1. Pre-registered version vs. readable versions

| File | Description |
|---|---|
| [`MSC_script.R`](https://github.com/Tamerabdalrazaq/MCS-Daniel-Psy-Ws/blob/main/MSC_script.R) | **The pre-registered version** of the analysis script, kept exactly as it was written for the pre-registration. |
| `MSC_script_<object>_<n>.R` (e.g. `MSC_script_c_n1.R`) | Re-organized versions of the same analysis, rewritten **for better readability and anonymization**. |

The analysis pipeline is the same, but there are some differences worth knowing about:

- **Single script vs. one script per condition.** The pre-registered script wraps the pipeline in a `run_script()` function and loops over both objects (`m`, `c`) and lags (`POL_EFFECT_ORDER` in `1:2`) in a single file. In the readable versions this is split into one script per study and lag (see below), and a third lag (n-3) was added.
- **Test subject.** The pre-registered script was written on the data of a single subject (`test_subject = 1`) and therefore *keeps only* that subject (`subjectnum == test_subject`). In the readable versions the same setting is used to *remove* that subject (`subjectnum != test_subject`), so that the reported analyses are run on the remaining participants and not on the data the pre-registration was based on.
- **Model saving.** The readable versions save every fitted `brms` model as an `.rds` file (named by lag, object, measure and saccade number), which is what the visualization script reads.
- **Anonymization.** Personal file paths in the pre-registered script are replaced by generic paths (`PATH_TO_MICROWAVES`, `PATH_TO_CLOCKS`).

## 2. Script naming convention

For convenience, the analysis is split into separate scripts, one per study and per "n-back" lag:

- **Lag (`n1`, `n2`, `n3`)** – the previous-target (PoL) effect from the target *n-1*, *n-2* and *n-3* trials back.
  - `n1` → n-1
  - `n2` → n-2
  - `n3` → n-3
- **Object (`m` / `c`)** – the study the script runs on.
  - `m` → microwaves (Study 1)
  - `c` → clocks (Study 2)

So, for example, `MSC_script_c_n1.R` runs the n-1 analysis on the clocks study.

The same settings are exposed at the top of each script:

```r
settings <<- list(
  object = 'c',          # 'm' for microwaves, 'c' for clocks
  POL_EFFECT_ORDER = 1,  # n-i effect (1, 2 or 3)
  ...
)
```

Before running, set `PATH_TO_MICROWAVES` / `PATH_TO_CLOCKS` (near the top of the script) to the location of the CSV files on your machine.

## 3. Exploratory analyses that were dropped

In the pre-registration we stated that, as exploratory analyses, we would examine:

1. **Dwell time** as a function of PoL
2. **Time spent near the previous target**

We decided to drop these two analyses for two reasons:

- **Theoretical relevance and space.** These analyses were theoretically less relevant to our main question, and given the space constraints of the journal, they were not important enough to be included.
- **Too many degrees of freedom.** The way these analyses were described in the AsPredicted pre-registration file left too many degrees of freedom (e.g., in how the measures are defined and computed), so they were not constrained enough to allow for a meaningful analysis.

We therefore do not report them.

The code for both is still included for transparency: in the pre-registered `MSC_script.R` they appear as "3. Dwell time as a function of PoL" and "5. Time spent near the previous target", and in the readable scripts as sections 5 and 6 (the numbering shifts because the readable scripts have additional sections). These sections are **not** part of the final analysis.

The other exploratory analysis ("Priming of last fixation vs. priming previous target") is not affected by this note.

## 4. Visualization script

We added a script that was not part of the pre-registered code:

- **`Visualizations.R`** – loads the saved model objects (`.RDS`) produced by the analysis scripts and generates the figures (posterior distributions of the PoL effect per saccade, for distance and direction/cosine measures, for microwaves and clocks and for n-1, n-2 and n-3).

Set `Path` (location of the saved models) and `saving_path` (output folder for figures) at the top of the script before running.

## 5. Supporting files

- **`Subject_outliers_Detector_function.R`** – helper function (`Sub_outlier_Detector`) for removing outlier observations per subject, using MAD (default) or mean/SD criteria. It is sourced by the analysis scripts.

## 6. Dataset

**`MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv`** is the processed fixation-level dataset for the microwaves study (one row per fixation).

### Columns inherited from the eye-tracking experiment

| Column | Meaning |
|---|---|
| `subjectnum` | Participant ID |
| `TRIAL_INDEX` | Trial number |
| `catcue` | Target category (cue) |
| `condition` | Target `present` / `absent` |
| `expected` | Correct response |
| `LAST_BUTTON_PRESSED` | Participant's response |
| `im_h`, `im_w` | Image height / width |
| `searcharray`, `img_name` | Image file name of the search array |
| `TRIAL_FIXATION_TOTAL` | Total number of fixations in the trial |
| `CURRENT_FIX_INDEX` | Index of the fixation within the trial |
| `CURRENT_FIX_X`, `CURRENT_FIX_Y` | Fixation location |
| `CURRENT_FIX_DURATION` | Fixation duration (ms) |

### Columns we added: object locations

We added the **locations of the target objects in each image** (the last five columns):

| Column | Meaning |
|---|---|
| `OBJECT_X1`, `OBJECT_Y1` | Top-left corner of the object's bounding box |
| `OBJECT_X2`, `OBJECT_Y2` | Bottom-right corner of the object's bounding box |
| `confidence` | Detection confidence for the bounding box |

These come from running an object detector on the search images. They make it possible to compute where the previous and current targets were located relative to each fixation (e.g., distance and direction of saccades relative to the previous target), and to filter out poor detections (`filter_bad_detections`) and very large objects (`filter_curr_big_objects`, `filter_prev_big_objects`).

For **target-absent** trials there is no object in the image, so all four coordinate columns and `confidence` are `N\A`. The scripts convert these to `NA` on load.
