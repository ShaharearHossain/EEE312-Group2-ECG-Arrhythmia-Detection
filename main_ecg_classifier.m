%% ====================================================================
% MAIN_ECG_CLASSIFIER  Execution Driver for EEE 312 DSP ECG Arrhythmia System
%
% Complete Pipeline Driver Script:
%   1. Loads reference dataset (.mat files) from data/ directory.
%   2. Runs Stage 3 ROC analysis to statistically determine opt_SDNN_cutoff.
%   3. Prompts user for a target .mat ECG file (or processes sample file).
%   4. Executes Stage 1 Pan-Tompkins algorithm for QRS complex detection.
%   5. Executes Stage 2 Feature Extraction (HR, SDNN, QRS Duration).
%   6. Evaluates Stage 4 Decision Rule Logic.
%   7. Displays structured console output dashboard and ECG plot with R, Q, S peaks.
%
% ====================================================================

clc;
clear;
close all;

fprintf('====================================================\n');
fprintf('  EEE 312 DSP PROJECT: ECG ARRHYTHMIA DETECTION SYSTEM\n');
fprintf('====================================================\n\n');

%% --------------------------------------------------------------------
%% STEP 1 & 2: Dataset Loading & Dynamic ROC Optimization (Stage 3)
%% --------------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

normal_folder     = fullfile(script_dir, 'data', 'normal');
arrhythmia_folder = fullfile(script_dir, 'data', 'arrhythmia');

fprintf('[STAGE 3] Loading reference dataset for ROC SDNN optimization...\n');

norm_files = dir(fullfile(normal_folder, '*.mat'));
arr_files  = dir(fullfile(arrhythmia_folder, '*.mat'));

fs_ref = 360;

% Extract SDNN features from Normal reference signals
normal_sdnn_list = [];
for k = 1:length(norm_files)
    mat_path = fullfile(normal_folder, norm_files(k).name);
    raw_data = load(mat_path);
    fnames   = fieldnames(raw_data);
    sig      = raw_data.(fnames{1});
    
    [R_p, Q_p, S_p] = detect_qrs_pan_tompkins(sig, fs_ref);
    [~, sdnn_val, ~, ~] = extract_qrs_features(R_p, Q_p, S_p, fs_ref, length(sig));
    if sdnn_val > 0
        normal_sdnn_list = [normal_sdnn_list; sdnn_val]; %#ok<AGROW>
    end
end

% Extract SDNN features from Arrhythmia reference signals
arrhythmia_sdnn_list = [];
for k = 1:length(arr_files)
    mat_path = fullfile(arrhythmia_folder, arr_files(k).name);
    raw_data = load(mat_path);
    fnames   = fieldnames(raw_data);
    sig      = raw_data.(fnames{1});
    
    [R_p, Q_p, S_p] = detect_qrs_pan_tompkins(sig, fs_ref);
    [~, sdnn_val, ~, ~] = extract_qrs_features(R_p, Q_p, S_p, fs_ref, length(sig));
    if sdnn_val > 0
        arrhythmia_sdnn_list = [arrhythmia_sdnn_list; sdnn_val]; %#ok<AGROW>
    end
end

% Execute Stage 3 ROC Analysis purely using extracted dataset features
[opt_SDNN_cutoff, sens, spec, thresholds] = perform_roc_analysis(normal_sdnn_list, arrhythmia_sdnn_list);
fprintf('[STAGE 3 COMPLETE] Dynamic Optimal SDNN Cutoff derived from dataset: %.2f ms\n\n', opt_SDNN_cutoff);

%% --------------------------------------------------------------------
%% STEP 3: Target ECG File Selection
%% --------------------------------------------------------------------
target_file = '';

if usejava('desktop')
    [fname, ppath] = uigetfile('*.mat', 'Select ECG Signal (.mat file) for Classification', script_dir);
    if isequal(fname, 0)
        fprintf('[INFO] No file selected. Exiting execution.\n');
        return; % Terminates the script immediately with no further output or plots
    end
    target_file = fullfile(ppath, fname);
else
    error('Interactive file selection requires desktop MATLAB environment.');
end

% Extract clean file name
[~, fname_only, ext_only] = fileparts(target_file);
selected_data_name = [fname_only, ext_only];

fprintf('Processing Target Signal File: %s\n', selected_data_name);

%% --------------------------------------------------------------------
%% STEP 4: Processing via Stage 1 (Pan-Tompkins) & Stage 2 (Feature Extraction)
%% --------------------------------------------------------------------
raw_struct = load(target_file);
fields = fieldnames(raw_struct);
ecg_signal = raw_struct.(fields{1});
fs = 360; % Standard MIT-BIH sampling rate

% Stage 1: QRS Detection
[R_peaks, Q_peaks, S_peaks] = detect_qrs_pan_tompkins(ecg_signal, fs);
plot_qrs_stages(ecg_signal, fs, R_peaks);

% Stage 2: Feature Extraction
[hr, sdnn, qrs_duration, rr_intervals_ms] = extract_qrs_features(R_peaks, Q_peaks, S_peaks, fs, length(ecg_signal));

%% --------------------------------------------------------------------
%% STEP 5: Stage 4 Classification Logic Evaluation
%% --------------------------------------------------------------------
[final_diagnosis, criteria] = classify_ecg_segment(hr, sdnn, qrs_duration, opt_SDNN_cutoff);

%% --------------------------------------------------------------------
%% STEP 6: Structured Console Output Dashboard
%% --------------------------------------------------------------------
fprintf('\n====================================\n');
fprintf('ECG SIGNAL ANALYSIS RESULT\n');
fprintf('====================================\n');
fprintf('Record Name   : %s\n', selected_data_name);
fprintf('Heart Rate    : %.2f bpm\n', hr);
fprintf('SDNN          : %.2f ms\n', sdnn);
fprintf('QRS Duration  : %.2f ms\n', qrs_duration);
fprintf('------------------------------------\n');
fprintf('FINAL DIAGNOSIS: %s\n', final_diagnosis);
fprintf('====================================\n\n');

%% --------------------------------------------------------------------
%% STEP 7: Signal & Peak Visualization Plot
%% --------------------------------------------------------------------
t = (0:length(ecg_signal)-1) / fs;

figure('Name', sprintf('ECG Signal Delineation [%s] - Diagnosis: %s', selected_data_name, final_diagnosis), ...
       'Color', 'w', 'Position', [150, 150, 950, 450]);

plot(t, ecg_signal, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Raw ECG Signal');
hold on;

if ~isempty(R_peaks)
    plot(t(R_peaks), ecg_signal(R_peaks), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'R Peaks');
end
if ~isempty(Q_peaks)
    plot(t(Q_peaks), ecg_signal(Q_peaks), 'bv', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'DisplayName', 'Q Peaks');
end
if ~isempty(S_peaks)
    plot(t(S_peaks), ecg_signal(S_peaks), 'gs', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'S Peaks');
end

xlabel('Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (mV / raw units)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('ECG Signal: %s [HR: %.1f bpm | SDNN: %.1f ms | QRS: %.1f ms] -> %s', ...
      selected_data_name, hr, sdnn, qrs_duration, final_diagnosis), ...
      'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
grid on;
grid minor;
lgd = legend('Location', 'northeast');
lgd.FontSize = 10;
xlim([0 min(10, max(t))]);
hold off;