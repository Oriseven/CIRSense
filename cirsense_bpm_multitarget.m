clear

diary off
% diary("cirsense_bpm_multitarget.txt");

% Initialize default parameters and configuration
params = build_default_params();
config = build_default_config();

params
config

% Get list of files from the data directory
files = dir(config.data_path);
files = files(~[files.isdir]);
[~, order] = sort([files.datenum], 'ascend');
files = files(order);
matchingFiles = string({files.name}).';

% Initialize results structure to store all measurements and comparisons
results = repmat(struct( ...
    "file_name", "", ...
    "respiration_rate_gt", NaN, ...
    "cirsense_bpm", NaN, ...
    "preprocess_time_per_frame", NaN, ...
    "dylign_time", NaN, ...
    "cirsense_total_time", NaN), 1, numel(matchingFiles));

% Initialize signal processing matrices for CSI processing
[~, used_index, used_sc, ~, Atime, ~] = initialization(params.used_sc,params.L);

% Process each selected file
for result_idx = 1:numel(matchingFiles)
    FileName = char(matchingFiles(result_idx));
    fprintf("\n");
    disp(FileName);

    data = load(config.data_path + FileName);
    csidata = data.csidata;     % CSI measurements
    t_slice = data.t_slice;     % Time slice information
    gt = data.gt;               % Ground truth respiration data
    gt_t = data.gt_t;           % Ground truth timestamps
    %%
    % Preprocess CSI data using the Domino algorithm
    tstart = tic;
    [csi_data_processed, cir_data_processed, ~] = Domino(csidata(:,:,params.used_antennas), length(params.used_antennas), used_index, used_sc, Atime, params);
    preprocess_time = toc(tstart);
    preprocess_time_per_frame = preprocess_time / size(csi_data_processed, 2);
    disp("preprocess time: " + num2str(preprocess_time_per_frame));

    % Calculate window length for smoothing based on sampling rate
    window_length = get_window_length(params.fs, params.smoothing_window, size(cir_data_processed, 2));

    % Calculate ground truth respiration rate from reference data
    respiration_rate_gt = bpm_est(gt, gt_t, 1, params);
    %%
    tstart = tic;
    [estimated_dis, h3, dylign_time] = cirsense(params.n_target, cir_data_processed, csi_data_processed, used_index, window_length, Atime, params);
    fprintf("\nEstimated distance by CIRSense: %s m\n", num2str(estimated_dis));
    respiration_rate_cir = bpm_est(h3(params.sel_target,:), t_slice, window_length, params);
    cirsense_total_time = toc(tstart);
    disp("Estimated bpm by CIRSense: " + num2str(respiration_rate_cir) + ...
        "  GT: " + num2str(respiration_rate_gt) + ...
        " Difference = " + num2str(abs(respiration_rate_gt - respiration_rate_cir)) + ...
        " using " + num2str(dylign_time) + " and " + num2str(cirsense_total_time));
    %%
    results(result_idx) = struct( ...
        "file_name", string(FileName), ...
        "respiration_rate_gt", respiration_rate_gt, ...
        "cirsense_bpm", respiration_rate_cir, ...
        "preprocess_time_per_frame", preprocess_time_per_frame, ...
        "dylign_time", dylign_time, ...
        "cirsense_total_time", cirsense_total_time);

    clear csi_com_all csi_time_com_all
end

result_table = struct2table(results);
csv_filename = 'multitarget_bmp_estimation_results_less1m.csv';
writetable(result_table, csv_filename);
fprintf('Results saved to %s\n', csv_filename);

diary off

function params = build_default_params()
% Initialize default parameters for signal processing and respiration estimation
params = struct();
params.used_antennas = 1;           % Index of antenna to use
params.nrx = 1;                     % Number of receiver antennas
params.L = -20:49;                  % Tap range for CIR
params.centraltap = 21;             % Central tap index
params.used_sc = [-984:-771, -765:-515, -509:-12, 12:509, 515:765, 771:984];  % Subcarrier indices to use
params.c = 3e8;                     % Speed of light (m/s)
params.los_dis = 0.6;               % Line-of-sight distance (m)
params.fs = 200;                    % Sampling frequency (Hz)
params.smoothing_window = 0.3;      % Smoothing window size in seconds
params.n_target = 2;                % Number of targets to detect
params.sel_target = 2;              % Index of target to focus on
% Parameters for Farsense algorithm
params.respiration_range = [10, 37];                    % Respiration rate range in BPM
params.respiration_range_hz = params.respiration_range / 60;  % Convert BPM to Hz
params.n_fft = 8192;                % FFT size for spectral analysis
end

function config = build_default_config()
config = struct();
config.data_path = "data\multitarget\target2_breathe\less1m\";
end

function window_length = get_window_length(fs, smoothing_window, ntime)
window_length = max(floor(min(fs * smoothing_window, ntime / 20) / 2) * 2 + 1, 5);
window_length = min(window_length, ntime - (mod(ntime, 2) == 0));
window_length = max(window_length, 3);
% Ensure window length is odd
if mod(window_length, 2) == 0
    window_length = window_length - 1;
end
end

function [estimated_dis, h, tend] = cirsense(n_target, cir_data, csi_data, used_index, smooth_window, Atime, params)

% Apply smoothing to reduce noise
smoothed_data = movmean(cir_data, smooth_window, 2);
H = movmean(csi_data, smooth_window, 2).';
% Initialize output arrays
estimated_dis = zeros(1, n_target);
h = zeros(n_target, size(smoothed_data, 2));
freq_factor = 1j * 2 * pi / 2048;

for t = 1:n_target
    % Find tap with maximum variance (indicates dynamic target)
    va = var(smoothed_data, [], 2);
    va(params.centraltap) = 0;
    [~, ti] = sort(va, "descend");
    approx_tap = ti(1);

    % Create range of candidate delays centered around the approximate tap
    tau_candidates_tap = linspace(approx_tap - params.centraltap - 0.5, approx_tap - params.centraltap + 0.5, 200);
    tstart = tic;
    phase_shift_candidate = exp(freq_factor * used_index.' * tau_candidates_tap);
    variance = zeros(1, numel(tau_candidates_tap));

    % Coarse search: evaluate candidate delays at larger intervals
    coarse_idx = 1:10:numel(tau_candidates_tap);
    phase_shift_pages = reshape(phase_shift_candidate(:, coarse_idx), 1, size(phase_shift_candidate, 1), []);
    H_shifted = permute(H(1:10:end, :) .* phase_shift_pages, [2,1,3]);
    h0 = pagemtimes(Atime(params.centraltap, :), H_shifted);
    variance(coarse_idx) = squeeze(var(h0, [], 2));
    [~, coarse_best] = max(variance(coarse_idx));

    % Fine search: evaluate candidate delays at finer resolution around best coarse result
    fine_start = max((coarse_best - 2) * 10 + 2, 1);
    fine_end = min(coarse_best * 10, numel(tau_candidates_tap));
    fine_idx = fine_start:fine_end;
    phase_shift_pages = reshape(phase_shift_candidate(:, fine_idx), 1, size(phase_shift_candidate, 1), []);
    H_shifted = permute(H(1:10:end, :) .* phase_shift_pages, [2,1,3]);
    h0 = pagemtimes(Atime(params.centraltap, :), H_shifted);
    variance(fine_idx) = squeeze(var(h0, [], 2));
    [~, local_best] = max(variance(fine_idx));

    % Select best delay and apply corresponding phase shift
    best_idx = fine_idx(local_best);
    H_sel = H .* phase_shift_candidate(:, best_idx).';
    tau_est = tau_candidates_tap(best_idx);

    % Extract breathing signal and calculate distance
    h(t, :) = mean(H_sel, 2).';
    estimated_dis(t) = tau_est/160e6*3e8;

    % Update H for next target
    H = H_sel;
    smoothed_data = Atime * H.';
end

tend = toc(tstart);
estimated_dis = (cumsum(estimated_dis) + params.los_dis) / 2;
[estimated_dis, order] = sort(estimated_dis);
h = h(order, :);
end

function respiration_rate = bpm_est(data, t_slice, smoothing_window, params)

[~, ntime] = size(data);
if smoothing_window > 1
    smoothdata = angle_projection_core(data, t_slice, ntime, params.respiration_range_hz, params.n_fft);
else
    smoothdata = data;
end

respiration_rate = estimate_respiration_rate(smoothdata, t_slice);
end

function [projected_signal, max_bnr] = angle_projection_core(H, t_slice, window_length, respiration_range_hz, n_fft)
% Calculate frequency axis based on signal duration
frequencies = build_frequency_axis(window_length, t_slice, n_fft);
% Find indices corresponding to respiration frequency range
respiration_indices = find(frequencies >= respiration_range_hz(1) & frequencies <= respiration_range_hz(2));
% Create range of projection angles
theta_range = linspace(0, pi, 100).';
% Project signal onto different angles
y = sin(theta_range) * real(H) + cos(theta_range) * imag(H);
y = normalize(y, 2);
% Calculate FFT for each projection
fft_length = max(n_fft, window_length);
signal_padded = [y, zeros(size(y, 1), fft_length - window_length)];
fft_result = fft(signal_padded, fft_length, 2);
% Calculate power in respiration frequency range
respiration_power = abs(fft_result(:, respiration_indices)).^2;
[~, local_peak_idx] = max(respiration_power, [], 2);
max_bin_index = respiration_indices(local_peak_idx);
total_energy = sum(abs(fft_result).^2, 2);
fft_result = reshape(fft_result.',[],1);
% Calculate energy at respiration frequency for each projection
respiration_energy = abs(fft_result(max_bin_index+[0:99]*fft_length)).^2;
% Calculate breathing-to-noise ratio (BNR) for each projection
BNRs = respiration_energy./total_energy;
% Select projection with highest BNR
[max_bnr, max_index] = max(BNRs);
projected_signal = y(max_index, :);
end

function frequencies = build_frequency_axis(window_length, t_slice, n_fft)
time_period = t_slice(end) - t_slice(1);
effective_fs = window_length / max(time_period, eps);
frequencies = (0:n_fft-1) * (effective_fs / n_fft);
end


function respiration_rate = estimate_respiration_rate(signal, t_slice)
signal = signal(:).';
max_lag = floor(numel(signal) / 3) * 2;
autocorr_results = calculate_autocorrelation(signal, max_lag);
[peaks, peak_locs] = findpeaks(autocorr_results, "MinPeakProminence", 0.1);

if isempty(peak_locs)
    respiration_rate = 0;
    return;
end

[~, max_peak_idx] = max(peaks);
peak_shift = peak_locs(max_peak_idx);
respiration_period = t_slice(peak_shift) - t_slice(1);
if respiration_period <= 0
    respiration_rate = 0;
else
    respiration_rate = 60 / respiration_period;
end
end

function autocorr_values = calculate_autocorrelation(signal, max_lag)
T = length(signal); % Length of the signal
y_i_mean = mean(signal); % Mean of the signal
% Initialize the autocorrelation values
autocorr_values = zeros(1, max_lag + 1);
% Calculate the denominator
denominator = sum((signal - y_i_mean).^2);
% Calculate the autocorrelation for each lag
for k = 0:max_lag
    numerator = sum((signal(k+1:T) - y_i_mean) .* (signal(1:T-k) - y_i_mean));
    autocorr_values(k+1) = numerator / denominator;
end
end


function [F, used_index, used_sc, FL, A, L] = initialization(used_sc, L)
F = ifftshift(ifftshift(dftmtx(2048), 1), 2);
% Select rows corresponding to used subcarriers
used_index = used_sc + 1024 + 1;
F = F(used_index, :);
% Select columns corresponding to time domain taps of interest
L = 1024 + L + 1;
used_index = used_sc + 1012 + 1;
FL = F(:, L);
% Create pseudo-inverse matrix for efficient CIR calculation
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
% Calculate initial CIR
cir_out = pagemtimes(Atime, csi_all);
% Find peak tap for each frame
[~, app_taps] = max(abs(cir_out));
Atime_center = Atime(params.centraltap, :);

for f = 1:size(csi_phase, 2)
    % Create range of candidate delays centered around peak tap
    tau_candidates_tap = linspace(app_taps(1, f, 1) - params.centraltap - 0.5, app_taps(1, f, 1) - params.centraltap + 0.5, 200);
    phase_shift_candidate = exp(freq_factor * used_sc.' * tau_candidates_tap);
    % Estimate sample timing offset
    [slop_new(f, 1), ~] = estimate_tau1_freq(csi_all(:, f, 1), tau_candidates_tap, phase_shift_candidate, Atime_center);

    for rx = 1:nrx
        % Apply phase correction to remove sample timing offset
        corrected_phase = csi_phase(:, f, rx) - freq_factor * 1j * slop_new(f, 1) * used_sc.';
        csi_out(:, f, rx) = csi_amp(:, f, rx) .* exp(1i .* corrected_phase);
        % Calculate corrected CIR
        cir_out(:, f, rx) = Atime * csi_out(:, f, rx);
        % Normalize to strongest path
        [~, max_index] = max(abs(cir_out(:, f, rx)));
        normalizer = cir_out(max_index, f, rx);
        csi_out(:, f, rx) = csi_out(:, f, rx) ./ normalizer;
        cir_out(:, f, rx) = cir_out(:, f, rx) ./ normalizer;
    end
end
end

function [tau1_est,i] = estimate_tau1_freq(H, tau_candidates_tap, phase_shift_candidate, Atime)
% Coarse search
H_shifted = H .* phase_shift_candidate(:,1:10:end);
[~,i] = max(abs(Atime*H_shifted));
% Fine search
fine_can = (i-2)*10+2:i*10;fine_can = fine_can(fine_can>0);
H_shifted = H .* phase_shift_candidate(:,fine_can);
[~,i] = max(abs(Atime*H_shifted));
tau1_est = tau_candidates_tap(fine_can(i));
i = fine_can(i);
end