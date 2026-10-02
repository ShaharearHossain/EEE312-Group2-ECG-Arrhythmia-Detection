function [hr, sdnn, qrs_duration, rr_intervals_ms] = extract_qrs_features(R_peaks, Q_peaks, S_peaks, fs, ~)
% EXTRACT_QRS_FEATURES Extract Quantitative Biomarkers from Detected QRS Complexes
%
% Syntax:
%    [hr, sdnn, qrs_duration, rr_intervals_ms] = extract_qrs_features(R_peaks, Q_peaks, S_peaks, fs)
%
% Inputs:
%    R_peaks  - Vector of detected R-peak sample indices
%    Q_peaks  - Vector of detected Q-peak sample indices
%    S_peaks  - Vector of detected S-peak sample indices
%    fs       - Sampling frequency in Hz (default: 360 Hz)
%
% Outputs:
%    hr              - Mean Heart Rate in beats per minute (bpm)
%    sdnn            - Standard deviation of RR intervals in milliseconds (ms)
%    qrs_duration    - Average QRS complex duration in milliseconds (ms)
%    rr_intervals_ms - Vector of consecutive RR intervals in milliseconds (ms)
%
% Formulae:
%    1. Heart Rate (HR):   bpm = 60 / mean(RR_sec)
%    2. SDNN (ms):         SDNN = std(RR_ms)
%    3. QRS Duration (ms): mean(|S_idx - Q_idx|) * (1000 / fs)

    % Default sampling frequency
    if nargin < 4 || isempty(fs)
        fs = 360;
    end

    % Ensure column vectors
    R_peaks = R_peaks(:);
    Q_peaks = Q_peaks(:);
    S_peaks = S_peaks(:);

    %% 1. RR Interval & Heart Rate (HR) Computation
    rr_sec = diff(R_peaks) / fs;
    rr_intervals_ms = rr_sec * 1000;
    hr = 60 / mean(rr_sec);
    sdnn = std(rr_intervals_ms);

    %% 2. Active QRS Complex Duration (ms)
    qrs_durations_samples = abs(S_peaks - Q_peaks);
    qrs_duration = mean(qrs_durations_samples) * (1000 / fs);

    %% 3. Rounding Output Features
    hr = round(hr, 2);
    sdnn = round(sdnn, 2);
    qrs_duration = round(qrs_duration, 2);
end