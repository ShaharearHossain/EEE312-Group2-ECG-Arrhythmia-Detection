function plot_qrs_stages(ecg_signal, fs, R_peaks)
% PLOT_QRS_STAGES  Visualize every Pan-Tompkins pipeline stage as subplots
%
% Syntax:
%    plot_qrs_stages(ecg_signal, fs, R_peaks)
%
% Inputs:
%    ecg_signal - Raw ECG signal vector (same signal passed to
%                 detect_qrs_pan_tompkins)
%    fs         - Sampling frequency in Hz
%    R_peaks    - (Optional) R-peak indices, e.g. from
%                 detect_qrs_pan_tompkins, to mark on the final subplot.
%                 Pass [] or omit if not available.

    if nargin < 2 || isempty(fs)
        fs = 360;
    end
    if nargin < 3
        R_peaks = [];
    end

    ecg_signal = ecg_signal(:);
    N = length(ecg_signal);

    %% ---------------- Stage 1: Normalize ----------------
    ecg_mean_removed = ecg_signal - mean(ecg_signal);
    max_amp = max(abs(ecg_mean_removed));
    if max_amp > 0
        ecg_norm = ecg_mean_removed / max_amp;
    else
        ecg_norm = ecg_mean_removed;
    end

    %% ---------------- Stage 2: Bandpass filter (5-15 Hz) ----------------
    lowcut = 5;
    highcut = 15;
    nyq = fs / 2;
    
    % 2nd-order Butterworth bandpass filter
    [b_bp, a_bp] = butter(2, [lowcut highcut] / nyq, 'bandpass');
    filtered_ecg = filtfilt(b_bp, a_bp, ecg_norm);

    %% ---------------- Stage 3: Derivative ----------------
    b_der = [-1 -2 0 2 1] / 8;
    derived_ecg = filter(b_der, 1, filtered_ecg);
    max_der = max(abs(derived_ecg));
    if max_der > 0
        derived_ecg = derived_ecg / max_der;
    end

    %% ---------------- Stage 4: Squaring ----------------
    squared_ecg = derived_ecg .^ 2;
    max_sq = max(abs(squared_ecg));
    if max_sq > 0
        squared_ecg = squared_ecg / max_sq;
    end

    %% ---------------- Stage 5: Moving window integration (~150 ms) ----------------
    win_samples = round(0.150 * fs);
    b_int = ones(1, win_samples) / win_samples;
    integrated_ecg = filter(b_int, 1, squared_ecg);
    max_int = max(abs(integrated_ecg));
    if max_int > 0
        integrated_ecg = integrated_ecg / max_int;
    end

    %% ---------------- Plot all stages as subplots ----------------
    t = (0:N-1) / fs;
    t_max = min(10, max(t));
    plot_mask = t <= t_max;

    figure('Name', 'Pan-Tompkins Pipeline - All Stages', 'Color', 'w', ...
           'Position', [100, 50, 950, 900]);

    subplot(3,2,1);
    plot(t(plot_mask), ecg_signal(plot_mask), 'k-', 'LineWidth', 1);
    title('Stage 0: Raw ECG Signal');
    ylabel('Amplitude'); grid on;

    subplot(3,2,2);
    plot(t(plot_mask), ecg_norm(plot_mask), 'b-', 'LineWidth', 1);
    title('Stage 1: Normalized (Mean Removed + Scaled)');
    ylabel('Amplitude'); grid on;

    subplot(3,2,3);
    plot(t(plot_mask), filtered_ecg(plot_mask), 'g-', 'LineWidth', 1);
    title('Stage 2: Bandpass Filtered (5-15 Hz)');
    ylabel('Amplitude'); grid on;

    subplot(3,2,4);
    plot(t(plot_mask), derived_ecg(plot_mask), 'm-', 'LineWidth', 1);
    title('Stage 3: Derivative Filtered');
    ylabel('Amplitude'); grid on;

    subplot(3,2,5);
    plot(t(plot_mask), squared_ecg(plot_mask), 'r-', 'LineWidth', 1);
    title('Stage 4: Squared');
    ylabel('Amplitude'); grid on;

    subplot(3,2,6);
    plot(t(plot_mask), integrated_ecg(plot_mask), 'c-', 'LineWidth', 1);
    hold on;
    if ~isempty(R_peaks)
        last_idx = find(plot_mask, 1, 'last');
        r_in_window = R_peaks(R_peaks <= last_idx);
        if ~isempty(r_in_window)
            plot(t(r_in_window), integrated_ecg(r_in_window), 'ko', ...
                 'MarkerFaceColor', 'y', 'MarkerSize', 6, 'DisplayName', 'R-peak');
        end
    end
    hold off;
    title('Stage 5: Moving-Window Integrated (~150 ms)');
    xlabel('Time (seconds)'); ylabel('Amplitude'); grid on;

    %% Backward-compatible global figure title
    title_text = 'Pan-Tompkins QRS Detection Pipeline - Stage by Stage';
    if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
        sgtitle(title_text, 'FontSize', 14, 'FontWeight', 'bold');
    else
        annotation('textbox', [0 0.94 1 0.05], 'String', title_text, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'FontSize', 12, 'FontWeight', 'bold');
    end
end