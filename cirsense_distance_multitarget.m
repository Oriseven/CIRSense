clear

diary off
% diary("cirsense_dis_2target.txt");

% Build default parameters and configuration
params = build_default_params();
config = build_default_config();

% Display parameters and config for verification
params
config

% Get list of data files in the specified directory, sorted by date
files = dir(config.data_path);
files = files(~[files.isdir]);
[~, order] = sort([files.datenum], 'ascend');
files = files(order);
matchingFiles = string({files.name}).';

% Initialize results structure to store estimation results for each file
results = repmat(struct( ...
    "file_name", "", ...
    "ground_truth_distance1", NaN, ...
    "ground_truth_distance2", NaN, ...
    "cirsense_dis1", NaN, ...
    "cirsense_dis2", NaN, ...
    "cirsense_error1", NaN, ...
    "cirsense_error2", NaN, ...
    "cirsense_time", NaN), 1, numel(matchingFiles));

% Initialize variables for signal processing (used subcarriers, etc.)
[~, used_index, used_sc, ~, Atime, ~] = initialization(params.used_sc, params.L);

%%
for result_idx = 1:numel(matchingFiles)
    FileName = char(matchingFiles(result_idx));
    fprintf("\n");
    disp(FileName);

    % Load CSI data, time slice, and ground truth distances from the file
    data = load(config.data_path + FileName);
    csidata = data.csidata;
    t_slice = data.t_slice;
    ground_truth_distance = data.ground_truth_distance;

    % Preprocess the CSI data using Domino algorithm
    tstart = tic;
    [csi_data_processed, cir_data_processed, ~] = Domino(csidata(:,:,params.used_antennas), length(params.used_antennas), used_index, used_sc, Atime, params);
    preprocess_time = toc(tstart);
    preprocess_time_per_frame = preprocess_time / size(csi_data_processed, 2);
    disp("Preprocessing time: " + num2str(preprocess_time_per_frame) + " s per frame");
    %%
    % Calculate window length for smoothing
    window_length = get_window_length(params.fs, params.smoothing_window, size(cir_data_processed, 2));

    % Perform distance estimation using CIRSense algorithm
    tstart = tic;
    csi_sel = movmean(csi_data_processed, window_length, 2);
    cir_sel = Atime * csi_sel;
    [cirsense_dis, var_value] = cirsense(params.n_target, cir_sel, csi_sel, used_index, Atime, window_length, params);
    cirsense_time = toc(tstart);

    % Calculate estimation errors
    cirsense_error1 = abs(cirsense_dis(1) - ground_truth_distance(1));
    cirsense_error2 = abs(cirsense_dis(2) - ground_truth_distance(2));

    % Display results
    fprintf("Ground truth distances: %.2f m, %.2f m\n", ground_truth_distance(1), ground_truth_distance(2));
    fprintf("CIRSense estimated distances: %.4f m, %.4f m\n", cirsense_dis(1), cirsense_dis(2));
    fprintf("CIRSense errors: %.4f m, %.4f m\n", cirsense_error1, cirsense_error2);
    fprintf("CIRSense processing time: %.4f s\n", cirsense_time);


    % Store results in the structure
    results(result_idx) = struct( ...
        "file_name", string(FileName), ...
        "ground_truth_distance1", ground_truth_distance(1), ...
        "ground_truth_distance2", ground_truth_distance(2), ...
        "cirsense_dis1", cirsense_dis(1), ...
        "cirsense_dis2", cirsense_dis(2), ...
        "cirsense_error1", cirsense_error1, ...
        "cirsense_error2", cirsense_error2, ...
        "cirsense_time", cirsense_time);

    % Clear large variables to free memory
    clear csi_data_processed cir_data_processed
end

% Convert results structure to table and save to CSV
result_table = struct2table(results);

csv_filename = 'multitarget_distance_estimation_results2.csv';
writetable(result_table, csv_filename);
fprintf('Results saved to %s\n', csv_filename);

diary off

function params = build_default_params()
% Build default parameters structure for the CIRSense algorithm
params = struct();
params.used_antennas = 1;  % Number of antennas used
params.nrx = 1;  % Number of receive antennas
params.L = -20:49;  % Tap range for CIR
params.centraltap = 21;  % Central tap index
params.used_sc = [-984:-771, -765:-515, -509:-12, 12:509, 515:765, 771:984];  % Used subcarriers
params.c = 3e8;  % Speed of light
params.los_dis = 0.6;  % Line-of-sight distance
params.fs = 500;  % Sampling frequency
params.smoothing_window = 1;  % Smoothing window in seconds
params.n_target = 2;  % Number of targets
params.sel_target = 1;  % Selected target index

% parameters for sigcan
params.threshold1 = 0.9;
params.threshold2 = 0.1;
params.N = 20;
params.subcarrier_spacing = 312.5e3;  % Subcarrier spacing
params.sc_idx = params.used_sc + 1012 +1;  % Subcarrier indices
params.sc_idx  = params.sc_idx(1:4:end);  % Downsampled indices
end

function config = build_default_config()
config = struct();
config.data_path = "data\multitarget\distance\";  % Path to data files
end

function window_length = get_window_length(fs, smoothing_window, ntime)
% Calculate the window length for moving average smoothing
% Ensures the window is odd, within bounds, and suitable for the data length
window_length = max(floor(min(fs * smoothing_window, ntime / 20) / 2) * 2 + 1, 5);
window_length = min(window_length, ntime - (mod(ntime, 2) == 0));
window_length = max(window_length, 3);
if mod(window_length, 2) == 0
    window_length = window_length - 1;
end
end

function [estimated_dis, h, tend] = cirsense(n_target, cir_data, csi_data, used_index, Atime, params)
% CIRSense algorithm for multitarget distance estimation
% Estimates distances to multiple targets using CIR and CSI data
% Inputs:
%   n_target: Number of targets
%   cir_data: Channel impulse response data
%   csi_data: Channel state information data
%   used_index: Indices of used subcarriers
%   Atime: Time-domain transformation matrix
%   params: Parameter structure (including window_length)
% Outputs:
%   estimated_dis: Estimated distances
%   h: Estimated channel responses
%   tend: Processing time

H = csi_data.';  % Transpose CSI data for processing
estimated_dis = zeros(1, n_target);  % Initialize distance estimates
h = zeros(n_target, size(cir_data, 2));  % Initialize channel responses
freq_factor = 1j * 2 * pi / 2048;  % Frequency factor for phase shifts
tstart = tic;  % Start timing

% Loop over each target
for t = 1:n_target
    % Compute variance across time for each tap
    va = var(cir_data, [], 2);
    va(params.centraltap) = 0;  % Zero out central tap
    [~, ti] = sort(va, "descend");  % Sort variances descending
    approx_tap = ti(1);  % Approximate tap with highest variance

    % Generate tau candidates around the approximate tap
    tau_candidates_tap = linspace(approx_tap - params.centraltap - 0.5, approx_tap - params.centraltap + 0.5, 200);
    phase_shift_candidate = exp(freq_factor * used_index.' * tau_candidates_tap);  % Phase shift candidates
    variance = zeros(1, numel(tau_candidates_tap));  % Variance array

    % Coarse search for best tau
    coarse_idx = 1:10:numel(tau_candidates_tap);
    phase_shift_pages = reshape(phase_shift_candidate(:, coarse_idx), 1, size(phase_shift_candidate, 1), []);
    H_shifted = permute(H(1:10:end, :) .* phase_shift_pages, [2,1,3]);
    h0 = pagemtimes(Atime(params.centraltap, :), H_shifted);
    variance(coarse_idx) = squeeze(var(h0, [], 2));
    [~, coarse_best] = max(variance(coarse_idx));  % Best coarse index

    % Fine search around the coarse best
    fine_start = max((coarse_best - 2) * 10 + 2, 1);
    fine_end = min(coarse_best * 10, numel(tau_candidates_tap));
    fine_idx = fine_start:fine_end;
    phase_shift_pages = reshape(phase_shift_candidate(:, fine_idx), 1, size(phase_shift_candidate, 1), []);
    H_shifted = permute(H(1:10:end, :) .* phase_shift_pages, [2,1,3]);
    h0 = pagemtimes(Atime(params.centraltap, :), H_shifted);
    variance(fine_idx) = squeeze(var(h0, [], 2));
    [~, local_best] = max(variance(fine_idx));  % Best fine index

    best_idx = fine_idx(local_best);  % Overall best index
    H_sel = H .* phase_shift_candidate(:, best_idx).';  % Apply best phase shift
    tau_est = tau_candidates_tap(best_idx);  % Estimated tau

    h(t, :) = mean(H_sel, 2).';  % Store channel response
    estimated_dis(t) = tau_est/160e6*3e8;  % Convert tau to distance

    % Update data for next target
    H = H_sel;
    cir_data = Atime * H.';
end

tend = toc(tstart);  % End timing
% Adjust distances: cumulative sum, add LOS, average, sort
estimated_dis = (cumsum(estimated_dis) + params.los_dis) / 2;
[estimated_dis, order] = sort(estimated_dis);
h = h(order, :);  % Sort channel responses accordingly
end

function [F, used_index, used_sc, FL, A, L] = initialization(used_sc, L)
% Initialize matrices and indices for signal processing
% Inputs:
%   used_sc: Used subcarrier indices
%   L: Tap indices
% Outputs:
%   F: DFT matrix subset
%   used_index: Adjusted subcarrier indices
%   used_sc: Used subcarriers
%   FL: F matrix subset for taps
%   A: Pseudoinverse for transformation
%   L: Adjusted tap indices

F = ifftshift(ifftshift(dftmtx(2048), 1), 2);  % DFT matrix with shifts
used_index = used_sc + 1024 + 1;  % Adjust indices
F = F(used_index, :);  % Subset of F
L = 1024 + L + 1;  % Adjust tap indices
used_index = used_sc + 1012 + 1;  % Another adjustment

FL = F(:, L);
A = (FL' * FL) \ FL';
end

function [csi_out, cir_out, slop_new] = Domino(csi_all, nrx, used_index, used_sc, Atime, params)
% Domino preprocessing algorithm for CSI data
% Corrects phase offsets and normalizes CSI data
% Inputs:
%   csi_all: Raw CSI data
%   nrx: Number of receive antennas
%   used_index: Used subcarrier indices
%   used_sc: Used subcarriers
%   Atime: Time-domain matrix
%   params: Parameters
% Outputs:
%   csi_out: Processed CSI data
%   cir_out: Time-domain CSI
%   slop_new: Slope estimates

s_frame = size(csi_all, 1);  % Number of frames
csi_all = csi_all(:, used_index, :);  % Select used subcarriers
csi_all = permute(csi_all, [2,1,3]);  % Permute dimensions

csi_amp = abs(csi_all);  % Amplitude
csi_phase = unwrap(angle(csi_all), [], 1);  % Unwrapped phase

csi_out = zeros(length(used_index), s_frame, nrx);  % Output CSI
slop_new = zeros(s_frame, 2);  % Slope estimates
freq_factor = 1j * 2 * pi / 2048;  % Frequency factor
cir_out = pagemtimes(Atime, csi_all);  % Time-domain transformation
[~, app_taps] = max(abs(cir_out));  % Approximate taps
Atime_center = Atime(params.centraltap, :);  % Central tap row

% Loop over frames
for f = 1:size(csi_phase, 2)
    % Generate tau candidates
    tau_candidates_tap = linspace(app_taps(1, f, 1) - params.centraltap - 0.5, app_taps(1, f, 1) - params.centraltap + 0.5, 200);
    phase_shift_candidate = exp(freq_factor * used_sc.' * tau_candidates_tap);  % Phase shifts
    [slop_new(f, 1), ~] = estimate_tau1_freq(csi_all(:, f, 1), tau_candidates_tap, phase_shift_candidate, Atime_center);  % Estimate slope

    % Process each receive antenna
    for rx = 1:nrx
        corrected_phase = csi_phase(:, f, rx) - freq_factor * 1j * slop_new(f, 1) * used_sc.';  % Correct phase
        csi_out(:, f, rx) = csi_amp(:, f, rx) .* exp(1i .* corrected_phase);  % Reconstruct CSI
        cir_out(:, f, rx) = Atime * csi_out(:, f, rx);  % Time-domain

        [~, max_index] = max(abs(cir_out(:, f, rx)));  % Find max tap
        normalizer = cir_out(max_index, f, rx);  % Normalizer
        csi_out(:, f, rx) = csi_out(:, f, rx) ./ normalizer;  % Normalize
        cir_out(:, f, rx) = cir_out(:, f, rx) ./ normalizer;  % Normalize time-domain
    end
end
end

function [tau1_est, i] = estimate_tau1_freq(H, tau_candidates_tap, phase_shift_candidate, Atime)
% Estimate the first tau (delay) in frequency domain
% Performs coarse and fine search for the best phase shift
% Inputs:
%   H: CSI vector
%   tau_candidates_tap: Candidate tau values
%   phase_shift_candidate: Corresponding phase shifts
%   Atime: Time-domain row
% Outputs:
%   tau1_est: Estimated tau
%   i: Index of best tau

% Coarse search
H_shifted = H .* phase_shift_candidate(:, 1:10:end);
[~, i] = max(abs(Atime * H_shifted));  % Find max after transformation

% Fine search around the coarse max
fine_can = (i-2)*10+2:i*10;
fine_can = fine_can(fine_can > 0);  % Ensure positive indices
H_shifted = H .* phase_shift_candidate(:, fine_can);
[~, i] = max(abs(Atime * H_shifted));  % Fine max
tau1_est = tau_candidates_tap(fine_can(i));  % Estimated tau
i = fine_can(i);  % Index
end