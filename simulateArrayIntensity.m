function [k,detection_mask] = simulateArrayIntensity(reflectivityImg, N, q, S, B)
%% FOR REFLECTIVITY ESTIMATION:
% simulate number of photon detections at each pixel

if ~exist('reflectivityImg','var') || isempty(reflectivityImg)
    imgSize = [500,500]; % image size in pixels
    a = ones(imgSize); % ground truth reflectivity of scene
else
    imgSize = size(reflectivityImg);
    a = reflectivityImg;
end
% figure; imshow(a,[]); colorbar; title('Ground truth reflectivity image');

% fixed dwell time: n pulses
if ~exist('N','var') || isempty(N)
    N = 1000;
end
if ~exist('q','var') || isempty(q)
    q = 0.17 + randn(imgSize)*0.05^2; % quantum efficiency (matrix over all pixels)
end
if ~exist('S','var') || isempty(S)
    h = fspecial('gaussian',500,250); % Gaussian beam falloff
    h = h/max(h(:));
    S = 0.09 * h; % average number of backreflected signal photons in one pulse repetition period (matrix over all pixels)
end
if ~exist('B','var') || isempty(B)
    B = 0.09*q; % average number of background photons reaching detector in one pulse repetition period
end


P0 = exp(-(q.*a.*S + B)); % probability of no detection in one pulse repetition period

% uniform random variables on [0,1]
% N 'events' per pixel
% at each pixel, count number of values > P0 --> photon detections, k(x,y)
% disp('Generating rvString...')
rvString = rand(imgSize(1), imgSize(2), N); 

% disp('Finding detections...')
detection_mask = bsxfun(@gt,rvString,P0); % 1 means detection, 0 means no detection

% disp('Finding number of detections at each pixel...')
k = sum(detection_mask,3); % simulated number of detections, k
% figure; imshow(k,[]); colorbar; title('Simulated number of photon detections');
% 
% detectionMask = sort(detection_mask,3,'descend');
% k_max = max(k(:));
% detectionMask = detectionMask(:,:,1:k_max);
