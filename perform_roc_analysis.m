function [opt_SDNN_cutoff, sensitivity, specificity, thresholds] = perform_roc_analysis(normal_sdnn_list, arrhythmia_sdnn_list)
% PERFORM_ROC_ANALYSIS  Statistical Threshold Optimization via ROC & Youden's J Index
%
% Syntax:
%   [opt_SDNN_cutoff, sensitivity, specificity, thresholds] = ...
%       perform_roc_analysis(normal_sdnn_list, arrhythmia_sdnn_list)
%
% Inputs:
%   normal_sdnn_list     - Vector of SDNN values from Normal reference ECG segments
%   arrhythmia_sdnn_list - Vector of SDNN values from Arrhythmia reference ECG segments
%
% Outputs:
%   opt_SDNN_cutoff      - Optimal decision cutoff maximizing Youden's J statistic
%   sensitivity          - Vector of Sensitivity (TPR) values across candidate cutoffs
%   specificity          - Vector of Specificity (1 - FPR) values across candidate cutoffs
%   thresholds           - Vector of candidate threshold cutoffs evaluated
%
% Algorithm:
%   1. Threshold Sweep across observed SDNN feature space.
%   2. Calculate TPR (Sensitivity) and FPR (1 - Specificity) at each threshold.
%   3. Maximize Youden's J index: J = Sensitivity + Specificity - 1.
%   4. Plot ROC curve with optimal operating point callout box.

    normal_sdnn_list = normal_sdnn_list(:);
    arrhythmia_sdnn_list = arrhythmia_sdnn_list(:);

    % Clean invalid entries
    normal_sdnn_list = normal_sdnn_list(~isnan(normal_sdnn_list) & normal_sdnn_list > 0);
    arrhythmia_sdnn_list = arrhythmia_sdnn_list(~isnan(arrhythmia_sdnn_list) & arrhythmia_sdnn_list > 0);

    P = length(arrhythmia_sdnn_list); % Positive class (Arrhythmia)
    N_norm = length(normal_sdnn_list);% Negative class (Normal)

    if P == 0 || N_norm == 0
        error('Input SDNN feature vectors must contain non-empty valid numerical data.');
    end

    %% 1. Threshold Sweep Setup
    all_features = [normal_sdnn_list; arrhythmia_sdnn_list];
    min_th = max(0, min(all_features) - 5);
    max_th = max(all_features) + 5;
    num_steps = 250;
    thresholds = linspace(min_th, max_th, num_steps);

    sensitivity = zeros(size(thresholds));
    specificity = zeros(size(thresholds));
    j_stats     = zeros(size(thresholds));

    %% 2. Metric Evaluation Across Cutoff Range
    for i = 1:length(thresholds)
        th = thresholds(i);

        % Arrhythmia is positive class (classified positive if SDNN >= th)
        TP = sum(arrhythmia_sdnn_list >= th);
        FN = sum(arrhythmia_sdnn_list < th);
        TN = sum(normal_sdnn_list < th);
        FP = sum(normal_sdnn_list >= th);

        sens = TP / (TP + FN);
        spec = TN / (TN + FP);

        sensitivity(i) = sens;
        specificity(i) = spec;
        j_stats(i)     = sens + spec - 1; % Youden's J statistic
    end

    %% 3. Youden's J Optimization for Optimal Cutoff
    [max_J, opt_idx] = max(j_stats);
    opt_SDNN_cutoff = thresholds(opt_idx);
    opt_sens = sensitivity(opt_idx);
    opt_spec = specificity(opt_idx);

    %% 4. Visualization: Plot ROC Curve
    one_minus_spec = 1 - specificity;

    figure('Name', 'ROC Curve & Youden J Optimal SDNN Cutoff Optimization', ...
           'Color', 'w', 'Position', [200, 200, 750, 550]);

    plot(one_minus_spec, sensitivity, 'b-', 'LineWidth', 2.5, 'DisplayName', 'ROC Curve');
    hold on;
    plot([0 1], [0 1], 'k--', 'LineWidth', 1.2, 'DisplayName', 'Random Chance Line');

    % Highlight Optimal Cutoff Point
    plot(1 - opt_spec, opt_sens, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
         'DisplayName', sprintf('Optimal Cutoff (%.2f ms)', opt_SDNN_cutoff));

    % Add Callout Text Box
    callout_str = sprintf('  Optimal Cutoff: %.2f ms\n  Sensitivity: %.2f%%\n  Specificity: %.2f%%\n  Youden J: %.4f', ...
                          opt_SDNN_cutoff, opt_sens * 100, opt_spec * 100, max_J);
    text(1 - opt_spec + 0.03, opt_sens - 0.05, callout_str, ...
         'BackgroundColor', [1 1 0.85], 'EdgeColor', 'r', 'FontSize', 10, 'FontWeight', 'bold');

    xlabel('1 - Specificity (False Positive Rate)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Sensitivity (True Positive Rate)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Receiver Operating Characteristic (ROC) Analysis for SDNN Cutoff', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    grid minor;
    lgd = legend('Location', 'southeast');
    lgd.FontSize = 10;
    axis([0 1 0 1]);
    hold off;

end