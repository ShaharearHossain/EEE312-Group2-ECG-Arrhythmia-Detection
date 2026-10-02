function [R_peaks, Q_peaks, S_peaks] = detect_qrs_pan_tompkins(ecg_signal, fs)
% DETECT_QRS_PAN_TOMPKINS  Pan-Tompkins QRS Complex Detector
%
% Syntax:
%   [R_peaks, Q_peaks, S_peaks] = detect_qrs_pan_tompkins(ecg_signal, fs)
%
% Inputs:
%   ecg_signal - Raw or pre-loaded 1D ECG signal vector (1xN or Nx1)
%   fs         - Sampling frequency in Hz (default: 360 Hz)
%
% Outputs:
%   R_peaks    - Indices of detected R-peaks
%   Q_peaks    - Indices of corresponding Q-peaks (preceding local minima)
%   S_peaks    - Indices of corresponding S-peaks (succeeding local minima)
%
% Algorithm Stages (Pan & Tompkins):
%   1. Preprocessing: Mean subtraction and peak amplitude normalization.
%   2. Bandpass Filtering (5-15 Hz): Attenuate baseline drift and T-wave/high-freq noise.
%   3. Derivative Operator: High-pass derivative to emphasize steep QRS slopes.
%   4. Squaring: Non-linear point transformation (y = x^2) emphasizing QRS energy.
%   5. Moving Window Integration: Smooth energy pulses over ~150 ms window.
%   6. Adaptive Thresholding & Peak Search: Refractory-constrained R, Q, S delineation.
%


    if nargin < 2 || isempty(fs)
        fs = 360; % Default MIT-BIH sampling rate
    end

    % Ensure input is a column vector
    ecg_signal = ecg_signal(:);
    N = length(ecg_signal);

    %% --------------------------------------------------------------------
    %% STAGE 1.1: Preprocessing (Mean Removal & Amplitude Normalization)
    %% --------------------------------------------------------------------
    ecg_mean_removed = ecg_signal - mean(ecg_signal);
    max_amp = max(abs(ecg_mean_removed));
    if max_amp > 0
        ecg_norm = ecg_mean_removed / max_amp;
    else
        ecg_norm = ecg_mean_removed;
    end

    %% --------------------------------------------------------------------
    %% STAGE 1.2: Bandpass Filtering (5 - 15 Hz)
    %% --------------------------------------------------------------------
    lowcut = 5;
    highcut = 15;
    nyq = fs / 2;
    
    % 2nd-order Butterworth bandpass filter
    [b_bp, a_bp] = butter(2, [lowcut highcut] / nyq, 'bandpass');
    filtered_ecg = filtfilt(b_bp, a_bp, ecg_norm);
    %% --------------------------------------------------------------------
    %% STAGE 1.3: Derivative Filtering
    %% --------------------------------------------------------------------
    % 5-point derivative operator: H(z) = (1/8)*(-z^-2 - 2z^-1 + 2z^1 + z^2)
    b_der = [-1 -2 0 2 1] / 8;
    derived_ecg = filter(b_der, 1, filtered_ecg);

    max_der = max(abs(derived_ecg));
    if max_der > 0
        derived_ecg = derived_ecg / max_der;
    end

    %% --------------------------------------------------------------------
    %% STAGE 1.4: Non-Linear Squaring Transformation
    %% --------------------------------------------------------------------
    squared_ecg = derived_ecg .^ 2;
    max_sq = max(abs(squared_ecg));
    if max_sq > 0
        squared_ecg = squared_ecg / max_sq;
    end

    %% --------------------------------------------------------------------
    %% STAGE 1.5: Moving Window Integration (~150 ms)
    %% --------------------------------------------------------------------
    win_samples = round(0.150 * fs); % 150 ms integration window
    b_int = ones(1, win_samples) / win_samples;
    integrated_ecg = filter(b_int, 1, squared_ecg);

    max_int = max(abs(integrated_ecg));
    if max_int > 0
        integrated_ecg = integrated_ecg / max_int;
    end

    %% --------------------------------------------------------------------
    %% STAGE 1.6: Adaptive Thresholding & Peak Search
    %% --------------------------------------------------------------------
    % Threshold computation: Combination of mean and peak amplitude
    mean_val = mean(integrated_ecg);
    max_val = max(integrated_ecg);
    TH = mean_val + 0.30 * (max_val - mean_val);

    % Refractory period: 200 ms
    refractory_samples = round(0.200 * fs);

    % Locate peak candidates in integrated signal
    raw_peaks = [];
    i = 1;
    while i <= N
        if integrated_ecg(i) > TH
            % Local maximum search within window
            search_end = min(N, i + refractory_samples);
            [~, max_rel_idx] = max(integrated_ecg(i:search_end));
            peak_pos = i + max_rel_idx - 1;
            raw_peaks = [raw_peaks; peak_pos]; %#ok<AGROW>
            i = peak_pos + refractory_samples;
        else
            i = i + 1;
        end
    end

    % Precise R-peak, Q-peak, and S-peak localization in normalized ECG signal
    R_peaks = [];
    Q_peaks = [];
    S_peaks = [];

    r_search_win = round(0.080 * fs); % 80 ms window for R-peak alignment
    q_win = round(0.060 * fs);        % 60 ms backward window for Q-peak
    s_win = round(0.060 * fs);        % 60 ms forward window for S-peak

    for k = 1:length(raw_peaks)
        p = raw_peaks(k);

        % Precise R-peak (local max in ECG signal near integrated peak)
        win_start = max(1, p - r_search_win);
        win_end = min(N, p + r_search_win);
        [~, rel_r] = max(ecg_norm(win_start:win_end));
        r_idx = win_start + rel_r - 1;

        % Q-peak (preceding minimum before R-peak)
        q_start = max(1, r_idx - q_win);
        [~, rel_q] = min(ecg_norm(q_start:r_idx));
        q_idx = q_start + rel_q - 1;

        % S-peak (succeeding minimum after R-peak)
        s_end = min(N, r_idx + s_win);
        [~, rel_s] = min(ecg_norm(r_idx:s_end));
        s_idx = r_idx + rel_s - 1;

        R_peaks = [R_peaks; r_idx]; 
        Q_peaks = [Q_peaks; q_idx]; 
        S_peaks = [S_peaks; s_idx]; 
    end
end
