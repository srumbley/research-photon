% Generates simulated photon detections.
%
% Input Parameters:
% - reflectivity: ground truth reflectivity image (n1 x n2).
% - depth: ground truth depth image (n1 x n2) in meters. If detection times
% not needed (only doing reflectivity estimation), input depth as [].
% - N: number of laser pulses.
% - A: expected number of detected signal photons from unity-reflectivity 
% pixel in one pulse repetition period, i.e. q.*S, where q is detection
% efficiency and S is average number of backreflected signal photons from
% unity-reflectivity pixel in one pulse repetition period.
% - B: expected number of detected background photons in one pulse
% repetition period (q.*B + d from paper).
% - Tr: pulse repetition period in seconds.
% - time_res: time resolution of detector in seconds (scalar).
% - pulse_args: two options...
%    1. provide pulse shape as histogram: 
%       pulse_args = {'pulse_hist', pulse_hist, pulse_bins};
%       - pulse_hist: pulse shape as histogram.
%       - pulse_bins: time bins for pulse_hist.
%    2. provide parameters to simulate truncated Gaussian pulse shape:
%       pulse_args = {'pulse_gauss', pulse_sigma, pulse_trunc};
%       - pulse_sigma: RMS pulse width (standard deviation of Gaussian).
%       - pulse_trunc: minimum and maximum pulse times (bounds of truncated
%       Gaussian), as [pulse_min, pulse_max].
%
% Outputs:
% - k: number of detections at each pixel. (n1 x n2)
% - t: 3D array (n1 x n2 x k_max) with detection times as first elements at
% each pixel, rest of elements are NaN.
% - detections: sparse matrix (n1*n2 x N) containing detection times.
% - signal_mask: ones indicate signal detections.

function [k,t,detections,signal_mask] = simulate_detections(reflectivity, depth, N, A, B, Tr, time_res, pulse_args)
disp('Getting detection mask...');
tic;
imSize = size(reflectivity);
numPixels = imSize(1)*imSize(2);

% probability of no detection in one pulse repetition period
P0 = exp(-(A.*reflectivity + B));

% at each pixel, randomly sample k (number of detections) from binomial
% distribution, where probability of success is 1-P0.
expected_counts = round(sum(N*(A(:).*reflectivity(:) + B(:))));
detection_mask = logical(spalloc(numPixels, N, expected_counts));
for i=1:N
    detection_mask(:,i) = rand(numPixels,1) > P0(:);
end
toc;
k = full(sum(detection_mask, 2));
k_max = max(k);
k = reshape(k, imSize);

if ~isempty(depth)

    % at each pixel, sample k detection times: signal times sampled from pulse 
    % shape distribution, noise sampled from uniform distribution on [0, Tr].
    disp('Getting detection times...');
    c = 3e8; % speed of light, m/ps
    t_truth = 2*depth/c;

    if ~exist('Tr','var') || isempty(Tr)
        Tr = 100e-9; % pulse repetition time in picoseconds
    end

    % determine whether each detection is due to signal or noise
    detections_rand = sprand(detection_mask);
    noise2signal_ratio = B./(A.*reflectivity + B);
    tic;
    signal_mask = logical(bsxfun(@gt,detections_rand,noise2signal_ratio(:)));
    noise_mask = detection_mask - signal_mask;

    % generate signal detection times = expected time of arrival + randomness from pulse shape.
    if (length(pulse_args) ~= 3)
        error('Pulse args input incorrectly.');
    else
        switch lower(pulse_args{1})
            case 'pulse_hist' % pulse shape as histogram
                pulse_hist = pulse_args{2};
                pulse_bins = pulse_args{3};
                
%                 pulse = normpdf(time_bins,0,sigma); % placeholder for testing

                pulse_hist = pulse_hist/sum(pulse_hist(:)); % normalize
                F = cumsum(pulse_hist);
                
                hist_bins = 0 : sum(pulse_hist)/100 : sum(pulse_hist); 
                % want fine resolution to map to every index in
                % pulse_bins...but if step size is smaller than
                % pulse_hist(1), then Finv(1)=0 --> time_bins(0) --> index
                % error. how to handle this? if Finv(1)=0, set Finv(1)=1 ?
                
                Finv = cumsum(histc(F, hist_bins));
                X = randi(length(Finv),sum(signal_mask(:)),1);
                rand_pulse_times = 10e-12 * time_bins(Finv(X));
                
            case 'pulse_gauss' % Gaussian pulse shape
                pulse_sigma = pulse_args{2};
                pulse_trunc = pulse_args{3};
                
                rand_pulse_times = TruncatedGaussian(pulse_sigma, pulse_trunc, [sum(signal_mask(:)),1]);
                
            otherwise
                error('Unrecognized pulse type.');
        end
    end

    signal_times = spdiags(t_truth(:),0,numPixels,numPixels) * signal_mask;
    signal_times(signal_mask) = signal_times(signal_mask) + rand_pulse_times;

    % generate noise detection times
    noise_times = sprand(noise_mask)*Tr;

    detections = signal_times + noise_times;

    % convert to detector time bins
    detections = ceil(detections/time_res);

    t = sort(detections,2,'descend');
    t = full(t(:,1:k_max));
    t = reshape(t,imSize(1),imSize(2),k_max);
    t(t==0) = NaN;
    toc;
end

end
