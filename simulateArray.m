% simulates data from array of photon-counting photodetectors.
%
% Arguments:
% - reflectivityImg: reflectivity of each pixel
% - depthImg: depth in meters at each pixel
% - N: number of pulse repetitions
% - qS: average number of detected signal photons from unity-reflectivity
% pixel in one pulse
% repetition period 
% - B: average number of detected background photons reaching detector in one
% pulse repetition period (including dark counts)
% - pulse: laser pulse shape represented as a histogram of number of
% photons emitted within each time bin
% - timeBins: time axis of pulse shape, in picoseconds. e.g. -250:25:250
% - Tr: pulse repetition time, in picoseconds.
%
% Returns:
% - k: number of photon detections at each pixel
% - t: 3D matrix, where size of the third dimension is equal to the maximum
% number of detections at a single pixel. t(i,j,1:k(i,j)) contains the
% times of all k(i,j) photon detections at pixel (i,j). The remaining
% t(i,j,k(i,j)+1:end) are filled in with NaN. 
% - t: cell array, where each cell t{i,j} contains a vector of the times of
% all photon detections at pixel (i,j). Time of arrival is measured from
% the last emitted pulse, in units of 8 picoseconds.

function [k,t,signal_mask] = simulateArray(reflectivity, depth, N, q, S, B, pulse, timeBins, Tr)
[k,detection_mask] = simulateArrayIntensity(reflectivity,N,q,S,B);

detection_mask = sort(detection_mask,3,'descend');
k_max = max(k(:));
detection_mask = detection_mask(:,:,1:k_max);

%% FOR RANGE ESTIMATION:
% simulate time of arrival data

imgSize = size(k);

if ~exist('depth','var') || isempty(depth)
    depth = peaks(imgSize(1)); % depth image
    depth = depth/max(abs(depth(:)));
    depth = depth + 2;
end
% figure
% surf(depthImg, 'LineStyle','none')
c = 3e-4; % meters/picosecond
t_actual = 2*depth/c;
% figure; surf(t_actual, 'LineStyle','none'); title('Ground truth TOA');
% figure; imshow(t_actual,[]); colorbar; title('Ground truth TOA');

if ~exist('timeBins','var') || isempty(timeBins)
    timeBins = -250:25:250; % time bins, width of 25 picoseconds each
end
if ~exist('Tr','var') || isempty(Tr)
    Tr = 100*1000; % pulse repetition time in picoseconds
end

% a few ways to do this:
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
%     pulse = normpdf(timeBins,0,sigma); % placeholder
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
    signal_times = timeBins(Finv(X)) + t_actual;
end
    
% determine whether each detection is due to signal or background
disp('Designating signal & background pixels...')
signal_mask = rand(imgSize(1),imgSize(2),k_max);
signal_mask = bsxfun(@le,signal_mask,(reflectivity.*q.*S)./(reflectivity.*q.*S + B));

% generate background detection times
disp('Generating background detections...')
bg_times = rand(imgSize(1),imgSize(2),k_max)*Tr;

t = (signal_mask.*signal_times + (1-signal_mask).*bg_times)./8;
t(detection_mask==0) = NaN;

% t_avg = nanmean(t,3); 
% figure; surf(t_avg,'LineStyle','none'); title('Simulated TOA (averaged over detections at each pixel)');
% figure; imshow(t_avg,[]); colorbar; title('Simulated TOA (averaged over detections at each pixel)');

% if (S>0)
%     t_signal = t;
%     t_signal(signal_mask==0) = NaN; % only signal detections
%     t_signal_avg = nanmean(t_signal,3); % average of signal detections
%     figure; imshow(t_signal_avg,[]); colorbar; title('Simulated signal-TOA (averaged over signal detections at each pixel)');
% end
end