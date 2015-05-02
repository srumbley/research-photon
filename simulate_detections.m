% generates simulated photon detections.
%
% Input Parameters:
% - reflectivity: ground truth reflectivity image (n1 x n2)
% - depth: ground truth depth image (n1 x n2)
% - N: number of laser pulses
% - A: number of detected signal photons from unity-reflectivity pixel in
% one pulse repetition period. equivalent to observation/sensing matrix for
% signal photons, i.e. q.*S, where q is detection efficiency and S is
% number of backreflected signal photons from unity-reflectivity pixel in
% one pulse repetition period.
% - B: number of detected background photons in one pulse repetition period
% - Tr: pulse repetition period in picoseconds
%
% Outputs:
% - detections: sparse matrix (n1*n2 x N) containing detection times.
% - k: number of detections at each pixel. (n1 x n2)
% - t: 3D array (n1 x n2 x k_max) with detection times as first elements at
% each pixel, rest of elements are NaN.

function [detections,k,t] = simulate_detections(reflectivity, depth, N, A, B, Tr, pulse, time_bins)
tic;
imgSize = size(reflectivity);

% probability of no detection in one pulse repetition period
P0 = exp(-(A.*reflectivity + B)); 

% uniform random variables on [0,1]
% N 'events' per pixel
% at each pixel, count number of values > P0 --> photon detections, k(x,y)
% disp('Generating rvString...')
rvString = rand(imgSize(1)*imgSize(2), N); 

% disp('Finding detections...')
detection_mask = bsxfun(@gt,rvString,P0(:)); % 1 means detection, 0 means no detection
detection_mask = sparse(detection_mask);

k = reshape(full(sum(detection_mask,2)), imgSize);
k_max = max(k(:));

% at pixel (i,j), first k(i,j) entries are 1, indicating detections.
% remaining entries from k(i,j)+1 to k_max are 0.
detection_mask_trunc = sort(detection_mask,2,'descend');
detection_mask_trunc = detection_mask_trunc(:,1:k_max);

toc;
tic;

c = 3e-4; % meters/picosecond
t_actual = 2*depth/c;

if ~exist('Tr','var') || isempty(Tr)
    Tr = 100*1000; % pulse repetition time in picoseconds
end

% generate signal detection times. a few ways to do this:
% 1. assume Gaussian curve: generate Gaussian random variables
if ~exist('pulse','var') || isempty(pulse)
    sigma = 100; % simulate Gaussian pulse shape

    % generate signal detection times
    disp('Generating signal detections...')
    gauss_signal = TruncatedGaussian(sigma,[-250,250],[imgSize(1),imgSize(2),k_max]);
    signal_times = bsxfun(@plus,gauss_signal,t_actual);

% 2. pulse shape as function: evaluate function over [-Tp/2, Tp/2], do
% method 3

% 3. pulse shape as histogram of photon counts: use cumsum(s)
else

    if ~exist('time_bins','var') || isempty(time_bins)
        time_bins = -250:25:250; % time bins, width of 25 picoseconds each
    end

%     pulse = normpdf(time_bins,0,sigma); % placeholder
    pulse = pulse/sum(pulse); % normalize
    F = cumsum(pulse); % cdf

    hist_bins = 0:sum(pulse)/100:sum(pulse); % want fine resolution to map to every index in time_bins.
                          % but if step size is smaller than s(1) -->
                          % Finv(1)=0 --> time_bins(0) --> index error.
                          % how to handle this?

    Finv = cumsum(histc(F,hist_bins));

    % generate signal detection times
    disp('Generating signal detections...')
    X = randi(length(Finv),imgSize(1),imgSize(2),k_max);
    signal_times = time_bins(Finv(X)) + t_actual;
end
    
% determine whether each detection is due to signal or background
disp('Designating signal & background pixels...')
signal_mask = rand(imgSize(1),imgSize(2),k_max);
signal_mask = bsxfun(@le,signal_mask,(A.*reflectivity)./(A.*reflectivity + B));

% generate background detection times
disp('Generating background detections...')
bg_times = rand(imgSize(1),imgSize(2),k_max)*Tr;

disp('Assigning detection times...')
times = (signal_mask.*signal_times + (1-signal_mask).*bg_times)./8;

t = times;
t(detection_mask_trunc==0) = NaN;

% find indices of detections in large matrix (n1 x n2 x N)
detection_inds = find(detection_mask);
[i,j] = ind2sub(size(detection_mask),detection_inds);

% create sparse matrix for detections
times = times(detection_mask_trunc==1);
detections = sparse(i,j,times,imgSize(1)*imgSize(2),N);
toc;
end
