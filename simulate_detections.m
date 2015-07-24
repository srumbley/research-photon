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

function [detections,signal_mask] = simulate_detections(reflectivity, depth, N, A, B, Tr, time_res, pulse_args)
% disp('Getting detection mask...');
% tic;
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
% toc;
% k = full(sum(detection_mask, 2));
% k_max = max(k);
% k = reshape(k, imSize);

if ~isempty(depth)

    % at each pixel, sample k detection times: signal times sampled from pulse 
    % shape distribution, noise sampled from uniform distribution on [0, Tr].
%     disp('Getting detection times...');
    c = 3e8; % speed of light, m/ps
    t_truth = 2*depth/c;

    if ~exist('Tr','var') || isempty(Tr)
        Tr = 100e-9; % pulse repetition time in picoseconds
    end

    % determine whether each detection is due to signal or noise
    detections_rand = sprand(detection_mask);
    noise2signal_ratio = B./(A.*reflectivity + B);
%     tic;
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
                step = 1e-5;
                if (pulse_hist(1) >= step)
                    pulse_hist(2:end+1) = pulse_hist;
                    pulse_hist(1) = 0;
                    pulse_bins(2:end+1) = pulse_bins;
                    pulse_bins(1) = pulse_bins(2);
                end
                pulse_hist = pulse_hist/sum(pulse_hist(:)); % normalize
                F = cumsum(pulse_hist);
                
                hist_bins = 0 : step : 1; 
                % want fine resolution to map to every index in
                % pulse_bins...but if step size is smaller than
                % pulse_hist(1), then Finv(1)=0 --> pulse_bins(0) --> index
                % error. how to handle this? if Finv(1)=0, set Finv(1)=1 ?
                
                Finv = cumsum(histc(F, hist_bins));
                X = randi(length(Finv)-2,sum(signal_mask(:)),1);
                rand_pulse_times = pulse_bins(Finv(X));
                
            case 'pulse_gauss' % Gaussian pulse shape
                pulse_sigma = pulse_args{2};
                pulse_trunc = pulse_args{3};
                
                rand_pulse_times = TruncatedGaussian(pulse_sigma, pulse_trunc, [sum(signal_mask(:)),1]);
                
            otherwise
                error('Unrecognized pulse type.');
        end
    end

    signal_times = spdiags(t_truth(:),0,numPixels,numPixels) * signal_mask;
    signal_times(signal_mask) = signal_times(signal_mask) + rand_pulse_times(:);

    % generate noise detection times
    noise_times = sprand(noise_mask)*Tr;

    detections = signal_times + noise_times;

    % convert to detector time bins
    detections = ceil(detections/time_res);

%     t = sort(detections,2,'descend');
%     t = full(t(:,1:k_max));
%     t = reshape(t,imSize(1),imSize(2),k_max);
%     t(t==0) = NaN;
%     toc;
end

% simulate crosstalk
% if (~isempty(p_ct))
%     ct_oob = 0;
%     ct_blocked = 0;
%     max_timebin = ceil(Tr/time_res);
%     ct_detections = sparse(numPixels, N);
%     ct_range = size(p_ct,3);
%     p_ct = reshape(p_ct, (2*imSize(1)-1)*(2*imSize(2)-1), ct_range);
%     
% %     tic;
%     for pulse = 1:N
%         frame = detections(:,pulse);
%         frame_im = reshape(frame, imSize);
%         [i,j,t] = find(frame_im);
%         ct_frame = nan(imSize);
% %         fprintf('pulse %d: %d detections, ', pulse, nnz(frame));
% %         tic;
%         for d = 1:nnz(frame)
%             is_ct = (rand((2*imSize(1)-1)*(2*imSize(2)-1), ct_range) < p_ct);
%             if sum(is_ct(:)) > 0
%                 time_offsets = -diag(ct_range:-1:1); % values reversed for taking min
%                 ct_offset = transpose(time_offsets * is_ct');
%                 ct_offset = min(ct_offset, [], 2);
%                 ct_offset = reshape(ct_offset, 2*imSize(1)-1, 2*imSize(2)-1);
% 
%                 [ct_di, ct_dj, ct_dt] = find(ct_offset);
%                 ct_di = ct_di - imSize(1);
%                 ct_dj = ct_dj - imSize(2);
%                 ct_dt = ct_dt + ct_range;
% 
%                 ct_frame_d = nan(imSize);
%                 ct_i = i(d) + ct_di;
%                 ct_j = j(d) + ct_dj;
% 
%                 oob_inds = [];
%                 for ii = 1:length(ct_i)
%                     if ct_i(ii) < 1 || ct_j(ii) < 1 || ct_i(ii) > imSize(1) || ct_j(ii) > imSize(2)
%                         oob_inds(end+1) = ii;
%                     end
%                 end
%                 if ~isempty(oob_inds)
%                     ct_oob = ct_oob + length(oob_inds);
%                     ct_i(oob_inds) = [];
%                     ct_j(oob_inds) = [];
%                     ct_dt(oob_inds) = [];
% 
%                     ct_indices = sub2ind(imSize,ct_i,ct_j);
%                     ct_frame_d(ct_indices) = mod(t(d) + ct_dt - 1, max_timebin-1) + 1;
%                     ct_blocked = ct_blocked + sum(~isnan(ct_frame(:))) + length(ct_indices);
%                     ct_frame = nanmin(ct_frame, ct_frame_d);
%                     ct_blocked = ct_blocked - sum(~isnan(ct_frame(:)));
%                 end
%             end
%         end
% %         toc;
% %         fprintf('%d crosstalk\n', sum(~isnan(ct_frame(:))));
%         ct_frame(isnan(ct_frame)) = 0;
%         ct_detections(:,pulse) = ct_frame(:);
%     end
% %     toc;
%     
%     detections = spfun(@(x) x - (max_timebin + 1), detections);
%     ct_detections = spfun(@(x) x - (max_timebin + 1), ct_detections);
%     detections = min(detections, ct_detections);
%     detections = spfun(@(x) x + max_timebin + 1, detections);
% end

end
