# ECG-Based Arrhythmia Detection System (MATLAB)

**EEE 312 – DSP Laboratory · Group 2** | Level 3, Term 1 | Dept. of EEE, BUET

| Student ID  |     Name                   |
|-------------|----------------------------|
| 2206100     | Joydip Das                 |
| 2206103     | Md Rakib                   |
| 2206105     | Md. Maheen Mosharrof Monon |
| 2206130     | Shaharear Hossain          |


## Overview
A MATLAB system that classifies a 10-second ECG segment as **NORMAL** or **ARRHYTHMIA** using QRS detection, three extracted features, and a ROC-optimized decision rule.

## Pipeline
1. **QRS detection** (`detect_qrs_pan_tompkins.m`): Pan–Tompkins algorithm (5–15 Hz band-pass, derivative, squaring, 150 ms moving-window integration, adaptive threshold) to locate R, Q and S peaks.
2. **Feature extraction** (`extract_qrs_features.m`): heart rate (bpm), SDNN (ms), QRS duration (ms).
3. **ROC analysis** (`perform_roc_analysis.m`): picks the SDNN cutoff that maximizes Youden's J (Sensitivity + Specificity − 1).
4. **Classification** (`classify_ecg_segment.m`): NORMAL only if **60 ≤ HR ≤ 100 bpm**, **SDNN < cutoff**, and **QRS ≤ 120 ms**; otherwise ARRHYTHMIA.

## Repository structure
```
main_ecg_classifier.m       # Run this
detect_qrs_pan_tompkins.m   # Stage 1
extract_qrs_features.m      # Stage 2
perform_roc_analysis.m      # Stage 3
classify_ecg_segment.m      # Stage 4
plot_qrs_stages.m           # Plots each pipeline stage
data/                       # Reference set: normal/ (48), arrhythmia/ (50)
test_data/                  # Demo set: test_normal/ (20), test_arrhythmia/ (20)
```

## How to run
**Requires:** MATLAB (desktop) with the Signal Processing Toolbox.

1. Clone or download this repo and set it as MATLAB's Current Folder.
2. Run `main_ecg_classifier` in the Command Window.
3. When prompted, select a `.mat` file from `test_data/`.

The program shows the ROC curve, the pipeline stages, and the ECG with detected R/Q/S peaks, and prints HR, SDNN, QRS duration and the final diagnosis.

