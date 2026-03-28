clear

diary off
% diary("cirsense_dis.txt");

params = build_default_params();
config = build_default_config();

params
config

files = dir(config.data_path);
files = files(~[files.isdir]);
[~, order] = sort([files.datenum], 'ascend');
files = files(order);

dis_str = ["10m","15m","20m","2m","2_25m","_3m","3_175m","4m","_5m","6m","7m"];
dis_val = [10, 15, 20, 2, 2.25, 3, 3.175, 4, 5, 6, 7];

matchingFiles = collect_matching_files(files, dis_str(config.distance_indices));
if isempty(matchingFiles)
    error("No files found matching the configured conditions.");
end

results = repmat(struct( ...
    "file_name", "", ...
    "ground_truth_distance", NaN, ...
    "cirsense_dis", NaN, ...
    "sigcan_distance", NaN, ...
    "cirsense_error", NaN, ...
    "sigcan_error", NaN, ...
    "cirsense_time", NaN, ...
    "sigcan_time", NaN), 1, numel(matchingFiles));

[~, used_index, used_sc, ~, Atime, ~] = initialization(params.used_sc, params.L);

%%
for result_idx = 1:numel(matchingFiles)
    FileName = char(matchingFiles(result_idx));
    fprintf("\n");
    disp(FileName);

    data = load(config.data_path + FileName);
    csidata = data.csidata;
    t_slice = data.t_slice;
    ground_truth_distance = data.ground_truth_distance;

    tstart = tic;
    [csi_data_processed, cir_data_processed, ~] = Domino(csidata(:,:,params.used_antennas), length(params.used_antennas), used_index, used_sc, Atime, params);
    preprocess_time = toc(tstart);
    preprocess_time_per_frame = preprocess_time / size(csi_data_processed, 2);
    disp("Preprocessing time: " + num2str(preprocess_time_per_frame) + " s per frame");
    %% CIRSense
    window_length = get_window_length(params.fs, params.smoothing_window, size(cir_data_processed, 2));

    tstart = tic;
    csi_sel = movmean(csi_data_processed, window_length, 2);
    cir_sel = Atime * csi_sel;
    [cirsense_dis, var_value] = cirsense(params.n_target, cir_sel, csi_sel, used_index, Atime, params);
    cirsense_time = toc(tstart);
    cirsense_error = abs(cirsense_dis - ground_truth_distance);

    fprintf("Ground truth distance: %.2f m\n", ground_truth_distance);
    fprintf("CIRSense estimated distance: %.4f m (error: %.4f m), ", cirsense_dis, cirsense_error);
    fprintf("CIRSense processing time: %.4f s\n", cirsense_time);
    %% SigCan
    tstart = tic;
    csi_smooth = movmean(csi_data_processed, 4, 1);
    csi_smooth = csi_smooth(1:4:end, :);
    csi_processed(1,:,:) = movmean(csi_smooth, window_length, 2);
    [~, sigcan_distance] = SigCan_ToF_Estimation(csi_processed, params);
    sigcan_time = toc(tstart);
    sigcan_error = abs(sigcan_distance - ground_truth_distance);

    fprintf("SigCan estimated distance: %.4f m (error: %.4f m), ", sigcan_distance, sigcan_error);
    fprintf("SigCan processing time: %.4f s\n", sigcan_time);

    results(result_idx) = struct( ...
        "file_name", string(FileName), ...
        "ground_truth_distance", ground_truth_distance, ...
        "cirsense_dis", cirsense_dis, ...
        "sigcan_distance", sigcan_distance, ...
        "cirsense_error", cirsense_error, ...
        "sigcan_error", sigcan_error, ...
        "cirsense_time", cirsense_time, ...
        "sigcan_time", sigcan_time);

    clear csi_data_processed cir_data_processed csi_processed
end

result_table = struct2table(results);

csv_filename = 'distance_estimation_results.csv';
writetable(result_table, csv_filename);
fprintf('Results saved to %s\n', csv_filename);

diary off

function params = build_default_params()
params = struct();
params.used_antennas = 1;
params.nrx = 1;
params.L = -20:49;
params.centraltap = 21;
params.used_sc = [-984:-771, -765:-515, -509:-12, 12:509, 515:765, 771:984];
params.c = 3e8;
params.los_dis = 0.6;
params.fs = 500;
params.smoothing_window = 1;
params.n_target = 1;
params.sel_target = 1;

% parameters for sigcan
params.threshold1 = 0.9;
params.threshold2 = 0.1;
params.N = 20;
params.subcarrier_spacing = 312.5e3;
params.sc_idx = params.used_sc + 1012 +1;
params.sc_idx  = params.sc_idx(1:4:end);
end

function config = build_default_config()
config = struct();
config.data_path = "data\distance\";
config.distance_indices = 1:11;
end

function matchingFiles = collect_matching_files(files, selected_dis_str)
file_names = string({files.name}).';
mask = false(size(file_names));

for idx = 1:numel(selected_dis_str)
    mask = mask | (contains(file_names, selected_dis_str(idx)));
end

matchingFiles = file_names(mask);
end

function window_length = get_window_length(fs, smoothing_window, ntime)
window_length = max(floor(min(fs * smoothing_window, ntime / 20) / 2) * 2 + 1, 5);
window_length = min(window_length, ntime - (mod(ntime, 2) == 0));
window_length = max(window_length, 3);
if mod(window_length, 2) == 0
    window_length = window_length - 1;
end
end

function [estimated_dis, h, tend] = cirsense(n_target, cir_data, csi_data, used_index, Atime, params)
H = csi_data.';
estimated_dis = zeros(1, n_target);
h = zeros(n_target, size(cir_data, 2));
freq_factor = 1j * 2 * pi / 2048;
tstart = tic;

for t = 1:n_target
    va = var(cir_data, [], 2);
    va(params.centraltap) = 0;
    [~, ti] = sort(va, "descend");
    approx_tap = ti(1);

    tau_candidates_tap = linspace(approx_tap - params.centraltap - 0.5, approx_tap - params.centraltap + 0.5, 200);
    phase_shift_candidate = exp(freq_factor * used_index.' * tau_candidates_tap);
    variance = zeros(1, numel(tau_candidates_tap));

    coarse_idx = 1:10:numel(tau_candidates_tap);
    phase_shift_pages = reshape(phase_shift_candidate(:, coarse_idx), 1, size(phase_shift_candidate, 1), []);
    H_shifted = permute(H(1:10:end, :) .* phase_shift_pages, [2,1,3]);
    h0 = pagemtimes(Atime(params.centraltap, :), H_shifted);
    variance(coarse_idx) = squeeze(var(h0, [], 2));
    [~, coarse_best] = max(variance(coarse_idx));

    fine_start = max((coarse_best - 2) * 10 + 2, 1);
    fine_end = min(coarse_best * 10, numel(tau_candidates_tap));
    fine_idx = fine_start:fine_end;
    phase_shift_pages = reshape(phase_shift_candidate(:, fine_idx), 1, size(phase_shift_candidate, 1), []);
    H_shifted = permute(H(1:10:end, :) .* phase_shift_pages, [2,1,3]);
    h0 = pagemtimes(Atime(params.centraltap, :), H_shifted);
    variance(fine_idx) = squeeze(var(h0, [], 2));
    [~, local_best] = max(variance(fine_idx));

    best_idx = fine_idx(local_best);
    H_sel = H .* phase_shift_candidate(:, best_idx).';
    tau_est = tau_candidates_tap(best_idx);

    h(t, :) = mean(H_sel, 2).';
    estimated_dis(t) = tau_est/160e6*3e8;

    H = H_sel;
    cir_data = Atime * H.';
end

tend = toc(tstart);
estimated_dis = (cumsum(estimated_dis) + params.los_dis) / 2;
[estimated_dis, order] = sort(estimated_dis);
h = h(order, :);
end

function [tau_est, distance_est] = SigCan_ToF_Estimation(CSI_data, params)
% 1. Static Paths and Dynamic Interference Cancellation

[csi_diff_sel, selected_sc] = path_cancellation(CSI_data, params);
if isempty(csi_diff_sel)
    tau_est = 0;
    distance_est = 0;
else
    % 2. Noise Reduction and Phase Processing
    [clean_phase] = noise_reduce(csi_diff_sel);

    % 3. ToF Estimation
    tau_est = clean_phase ./ ([0:selected_sc(2)-selected_sc(1)].' * params.subcarrier_spacing);
    tau_est = mean(tau_est(tau_est > 0));
    distance_est = (tau_est * params.c + params.los_dis) / 2; % round-trip
end
end

function [csi_sel, selected_sc] = path_cancellation(CSI, params)
% Phase Linearity & Signal Similarity Metrics
[nAnt, nSc, nTime] = size(CSI);
CSI_diff = CSI(:, :, 2:end) - CSI(:, :, 1); % Time-domain difference for static path removal

sc_idx = params.sc_idx; % Center subcarrier indices
sc_20 = floor(nSc/8);
sc_range = 1:sc_20:nSc-sc_20;
time_sel = {};
metric2 = zeros(nSc, 2);
valid_ranges = {};
valid_performace = {};
num = 0;

for a = 1:nAnt
    for band = 1:length(sc_range(1:end-1))
        num = num + 1;
        % Metric 1: Phase Linearity
        ref_f = sc_range(band) + sc_20/2;
        sel_sc_range = sc_range(band):sc_range(band+1)-1;
        delta_f = (sc_idx(sel_sc_range) - sc_idx(ref_f)) .* params.subcarrier_spacing;
        delta_f(1+sc_20/2) = [];  % frequency interval between different subcarriers
        sc_phase = squeeze(angle(CSI_diff(a, sel_sc_range, :)));
        sc_phase = unwrap(sc_phase);
        phase_diff = sc_phase - sc_phase(1+sc_20/2, :);
        phase_diff(1+sc_20/2, :) = []; % phase difference of the time-domain difference vector

        corr_coeffs = zeros(nTime-1, 1);
        for t = 1:nTime-1
            [R, ~] = corrcoef(delta_f, phase_diff(:, t));
            corr_coeffs(t) = R(1, 2); % correlation coefficient -> linearity
        end
        time_sel{end+1} = find(abs(corr_coeffs) > params.threshold1);  % select time domain frame based on metric 1

        % Metric 2: Signal Similarity
        for s = sel_sc_range
            objective = @(params) mean(abs(CSI(a, ref_f, time_sel{num}+1) - CSI(a, s, time_sel{num}+1) .* exp(1j * 2 * pi * params(1))));
            options = optimset('Display', 'off');
            [~, D] = fminsearch(objective, 0, options);
            metric2(s, a) = D;
        end

        [valid_ranges{end+1}, valid_performace{end+1}] = findConsecutiveSubcarriers(sc_20, metric2(sel_sc_range, a), params.N, params.threshold2);
    end
end

best_performance = inf;
select_antenna = 0;
not_valid = 1;

for a = 1:nAnt
    for band = 1:length(sc_range(1:end-1))
        n = (a-1) * length(sc_range(1:end-1)) + band;
        sel_sc_range = sc_range(band):sc_range(band+1)-1;

        if isempty(valid_ranges{n})
            selected_sc = [1, nSc];
        else
            [performance, i] = min(valid_performace{n});
            if performance < best_performance % for metric 2, smaller is better
                select_antenna = a;
                sel_time = n;
                current_candidate = valid_ranges{n};
                selected_sc = sel_sc_range(current_candidate(i, :));
                best_performance = performance;
            end
            not_valid = 0;
        end
    end
end

if not_valid
    csi_sel = [];
else
    csi_sel = squeeze(CSI_diff(select_antenna, selected_sc(1):selected_sc(2), time_sel{sel_time}));
end
end

function [clean_phase] = noise_reduce(CSI_sel)
% Noise Reduction
[nSc, ~] = size(CSI_sel);
ref_f = 1;
clean_phase = zeros(nSc, 1);

for s = 1:nSc
    objective = @(params) sum(abs(CSI_sel(ref_f, :) - CSI_sel(s, :) .* exp(1j * 2 * pi * params(1))).^2);
    options = optimset('Display', 'off');
    [clean_phase(s), ~] = fminsearch(objective, 0, options);
end
end

function [valid_ranges, valid_performace] = findConsecutiveSubcarriers(subcarriers, metric2, N, threshold2)
valid_ranges = [];
valid_performace = [];

for i = 1:(subcarriers - N + 1)
    % Extract the metrics for the current window of size N
    metric2_window = metric2(i:i+N-1);

    % Check if all metrics in the window exceed their thresholds
    if all(metric2_window < threshold2)
        % Store the valid range
        valid_ranges = [valid_ranges; i, i + N - 1];
        valid_performace = [valid_performace; mean(metric2_window)];
        % Use the mean value to represent the performance of this candidate
    end
end
end

function [F, used_index, used_sc, FL, A, L] = initialization(used_sc, L)
F = ifftshift(ifftshift(dftmtx(2048), 1), 2);
used_index = used_sc + 1024 + 1;
F = F(used_index, :);
L = 1024 + L + 1;
used_index = used_sc + 1012 + 1;

FL = F(:, L);
A = (FL' * FL) \ FL';
end

function [csi_out, cir_out, slop_new] = Domino(csi_all, nrx, used_index, used_sc, Atime, params)
s_frame = size(csi_all, 1);
csi_all = csi_all(:, used_index, :);
csi_all = permute(csi_all, [2,1,3]);

csi_amp = abs(csi_all);
csi_phase = unwrap(angle(csi_all), [], 1);

csi_out = zeros(length(used_index), s_frame, nrx);
slop_new = zeros(s_frame, 2);
freq_factor = 1j * 2 * pi / 2048;
cir_out = pagemtimes(Atime, csi_all);
[~, app_taps] = max(abs(cir_out));
Atime_center = Atime(params.centraltap, :);

for f = 1:size(csi_phase, 2)
    tau_candidates_tap = linspace(app_taps(1, f, 1) - params.centraltap - 0.5, app_taps(1, f, 1) - params.centraltap + 0.5, 200);
    phase_shift_candidate = exp(freq_factor * used_sc.' * tau_candidates_tap);
    [slop_new(f, 1), ~] = estimate_tau1_freq(csi_all(:, f, 1), tau_candidates_tap, phase_shift_candidate, Atime_center);

    for rx = 1:nrx
        corrected_phase = csi_phase(:, f, rx) - freq_factor * 1j * slop_new(f, 1) * used_sc.';
        csi_out(:, f, rx) = csi_amp(:, f, rx) .* exp(1i .* corrected_phase);
        cir_out(:, f, rx) = Atime * csi_out(:, f, rx);

        [~, max_index] = max(abs(cir_out(:, f, rx)));
        normalizer = cir_out(max_index, f, rx);
        csi_out(:, f, rx) = csi_out(:, f, rx) ./ normalizer;
        cir_out(:, f, rx) = cir_out(:, f, rx) ./ normalizer;
    end
end
end

function [tau1_est, i] = estimate_tau1_freq(H, tau_candidates_tap, phase_shift_candidate, Atime)
H_shifted = H .* phase_shift_candidate(:, 1:10:end);
[~, i] = max(abs(Atime * H_shifted));
fine_can = (i-2)*10+2:i*10;
fine_can = fine_can(fine_can > 0);
H_shifted = H .* phase_shift_candidate(:, fine_can);
[~, i] = max(abs(Atime * H_shifted));
tau1_est = tau_candidates_tap(fine_can(i));
i = fine_can(i);
end