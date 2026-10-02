function [diagnosis, criteria_met] = classify_ecg_segment(hr, sdnn, qrs_duration, opt_SDNN_cutoff)
% CLASSIFY_ECG_SEGMENT  Stage 4 Decision Rule Logic for ECG Classification
%
% Syntax:
%    [diagnosis, criteria_met] = classify_ecg_segment(hr, sdnn, qrs_duration, opt_SDNN_cutoff)
%
% Inputs:
%    hr              - Heart Rate in beats per minute (bpm)
%    sdnn            - RR interval standard deviation in milliseconds (ms)
%    qrs_duration    - QRS duration in milliseconds (ms)
%    opt_SDNN_cutoff - Statistically optimized SDNN threshold derived via ROC
%
% Outputs:
%    diagnosis       - String: 'NORMAL' or 'ARRHYTHMIA'
%    criteria_met    - Struct containing boolean evaluation of individual rules:
%                        .hr_normal     (60 <= hr <= 100)
%                        .sdnn_normal   (sdnn < opt_SDNN_cutoff)
%                        .qrs_normal    (qrs_duration <= 120)
%
% Decision Logic:
%    NORMAL     : Requires ALL THREE conditions to hold true simultaneously.
%    ARRHYTHMIA : Triggered if ANY condition fails.


    % Evaluate individual boolean conditions using dataset-derived ROC cutoff
    cond_hr   = (hr >= 60.0) && (hr <= 100.0);
    cond_sdnn = (sdnn < opt_SDNN_cutoff);
    cond_qrs  = (qrs_duration <= 120.0);

    criteria_met = struct();
    criteria_met.hr_normal   = cond_hr;
    criteria_met.sdnn_normal = cond_sdnn;
    criteria_met.qrs_normal  = cond_qrs;

    % Final Classification: ALL three criteria must pass for NORMAL
    if cond_hr && cond_sdnn && cond_qrs
        diagnosis = 'NORMAL';
    else
        diagnosis = 'ARRHYTHMIA';
    end
end