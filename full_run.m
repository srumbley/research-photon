load('gt_data/tank2');
I_truth = imresize(I_truth,0.25);
D_truth = D_truth(1:4:end,1:4:end);

im_size = size(I_truth);
num_pixels = im_size(1)*im_size(2);
array_size = [32,32];
array_reps = ceil(im_size./array_size);

N = 200;
q_mean = 0.35;
d_mean = 0.0022;
h = fspecial('gaussian',array_size,0.9*array_size(1)); % Gaussian beam falloff
h = h/max(h(:));
S = 0.17 * h;
b = 0.009;

c = 3e8;
Tr = 1e-6;
time_res = 1e-9;
max_timebin = ceil(Tr/time_res);
pulse_func = @(x) ((x-80e-12).^4).*exp(-(x-80e-12)./(40e-12));
pulse_bins = 90e-12:1e-12:1e-9;
pulse_hist = pulse_func(pulse_bins);
pulse_hist = pulse_hist/sum(pulse_hist);
pulse_args = {'pulse_hist', pulse_hist, pulse_bins};

load('crosstalk_probability2');
% p_ct_sim = p_ct_sim/sum(p_ct_sim(:))*0.15/q_mean/0.4111;
p_ct = p_ct/sum(p_ct(:));
p_ct = p_ct*0.15/q_mean/sum(sum(sum(p_ct(:,16:48,16:48)))); % 16:48 only for array size 32x32

im2array_inds = reshape((1:num_pixels),im_size);
im2array_inds = im2col(im2array_inds,array_size,'distinct');

% SPIRAL parameters
I1_tau = 3.5; I1_ainit = 100; I1_maxiter = 100;
I_penalty = 'tv'; I_noisetype = 'binomial';
D1_tau = 1; D1_ainit = 1; D1_maxiter = 50; 
D_penalty = 'tv'; D_noisetype = 'gaussian';
I2_tau = 4; I2_ainit = 10; I2_maxiter = 100;
D2_tau = 0.7; D2_ainit = 1; D2_maxiter = 50; 

filedir = 'results/full/tank2_2/';
save([filedir, 'common_vars']);

%% Begin loop - run repeated trials
for iter=1:1
    %% Simulate detections
    
    % generate randomized array parameters (quantum efficiency, dark counts, crosstalk probability, dead pixels)
    q = q_mean + randn(array_size)*0.0167;
    d = d_mean + randn(array_size)*0.0004;
    
%     p_ct = bsxfun(@plus,p_ct_sim,randn(1,size(p_ct_sim,2),size(p_ct_sim,3))*0.0002^2);
%     p_ct(p_ct<0) = 0;
%     p_ct(:,array_size(2),array_size(1)) = 0;
    
    deadpix = rand(array_size)>0.001;
    
    signal_level = 0.5*q.*S; % 0.5*q because incoming photons are attenuated by bandpass filter
    noise_level = 0.5*q.*b + d;
    noise_level_im = repmat(noise_level,array_reps);
    signal_level_im = repmat(signal_level,array_reps);
    
    dpr = sparse(num_pixels, N);
    signal_mask = logical(sparse(num_pixels, N));
    detections = sparse(num_pixels, N);
    tic;
    for n = 1:array_reps(1)
        fprintf('scan row %d: ',n);
        i_array = (n-1)*array_size(1)+1 : min(n*array_size(1), im_size(1));
        for m = 1:array_reps(2)
            fprintf('%d ', m); 
            j_array = (m-1)*array_size(2)+1 : min(m*array_size(2), im_size(2));
            [dpr_array, signal_mask_array] = simulate_detections(I_truth(i_array,j_array), D_truth(i_array,j_array), N,signal_level,noise_level,Tr,time_res,pulse_args);
            dpr_array(~deadpix,:) = 0;
            signal_mask_array(~deadpix,:) = 0;
            detections_array = simulate_crosstalk(dpr_array, p_ct, ceil(Tr/time_res), q, deadpix);
            array_inds = im2array_inds(:, array_reps(1)*(m-1)+n);
            array_inds(array_inds==0) = [];
            dpr(array_inds,:) = dpr_array;
            signal_mask(array_inds,:) = signal_mask_array;
            detections(array_inds,:) = detections_array;
        end
        fprintf('\n');
    end
    toc;
    ct_mask = ((dpr==0) & (detections~=0)) | (dpr > detections);
    signal_mask(detections < dpr) = 0;
    noise_mask = (dpr ~= 0) & (~ct_mask) & (~signal_mask);

    [k,t] = get_kt(detections,im_size);

    fprintf(['fraction of pixels with no detections: ' num2str(sum(k(:)==0)/length(k(:))) '\n']);
    fprintf(['max # detections: ' num2str(max(k(:))) '\n']);
    fprintf(['mean # detections: ' num2str(mean(k(:))) '\n']);
    fprintf(['percentage crosstalk: ' num2str(nnz(ct_mask)/(nnz(signal_mask)+nnz(noise_mask))) '\n']);
    
    save([filedir, sprintf('%d_simulated_data',iter)],'q','d','p_ct','deadpix','detections','signal_mask','ct_mask','noise_mask','k','t');

    %% Compute crosstalk probability for each detection (based on previous detections)
    
    ct_prob = sparse(num_pixels, N);
    for n = 1:array_reps(1)
        i_array = (n-1)*array_size(1)+1 : min(n*array_size(1), im_size(1));
        for m = 1:array_reps(2)
            j_array = (m-1)*array_size(2)+1 : min(m*array_size(2), im_size(2));

            array_inds = im2array_inds(:, array_reps(1)*(m-1)+n);
            array_inds(array_inds==0) = [];

            ct_prob(array_inds,:) = compute_ct_prob(detections(array_inds,:), p_ct);
        end
    end
    ct_prob_im = reshape(full(sum(ct_prob,2)),im_size);
    
    %% Initial intensity estimation

    fprintf('\n Estimating intensity... \n');

    % max likelihood estimation
    I_ML = (log(N./(N-k))-noise_level_im-ct_prob_im)./signal_level_im;
    
    % SPIRAL parameters
    I_A = @(x) signal_level_im.*x; 
    I_AT = @(x) signal_level_im.*x; 

    I_MAP = SPIRALTAP_NEW(k,I_A,I1_tau,N, ...
        'noisetype',I_noisetype, ... % approx for binomial
        'penalty',I_penalty, ...
        'maxiter',I1_maxiter,... % `120
        'Initialization',I_ML./max(I_ML(:)),...
        'AT',I_AT,...
        'monotone',1,...
        'miniter',1,...
        'stopcriterion',3,...
        'tolerance',1e-8,...
        'alphainit',I1_ainit,...
        'alphaaccept',1e80,...
        'logepsilon',1e-10 + noise_level_im + ct_prob_im,...
        'saveobjective',1,...
        'savereconerror',1,...
        'savecputime',1,...
        'savesolutionpath',0,...
        'truth',I_truth,...
        'verbose',10);

    I_psnr = psnr(I_MAP,I_truth)
    figure;imshow(I_MAP)
    
    % compute offset between I_MAP and I_truth, and psnr when I_MAP adjusted for
    % offset.
    [counts,centers] = hist(I_MAP(:)-I_truth(:),100);
    [~,offset_ind] = max(counts);
    I_offset = centers(offset_ind);
    I_offset_psnr = psnr(I_MAP-I_offset,I_truth)
    
    %% Censor background noise (ROM filtering)
    
    rom = get_rom(t,2);
    threshold = 5*(2+size(p_ct,1))*(noise_level_im)./(I_MAP.*signal_level_im + noise_level_im + ct_prob_im);
    threshold(isnan(threshold)) = max_timebin;
    d_censored = censor(detections,rom,threshold);
    
    signal_censored = nnz((~d_censored)&signal_mask)/nnz(signal_mask);
    noise_censored = nnz((~d_censored)&noise_mask)/nnz(noise_mask);
    ct_censored = nnz((~d_censored)&ct_mask)/nnz(ct_mask);
    fprintf('accuracy = %f \n', nnz(((d_censored~=0)==signal_mask)&detections)/nnz(detections));
    fprintf('fraction of background censored = %f \n', nnz((d_censored==0)&noise_mask)/nnz(noise_mask));
    fprintf('fraction of crosstalk censored = %f \n', nnz((d_censored==0)&ct_mask)/nnz(ct_mask));
    fprintf('fraction of signal retained = %f \n', nnz(d_censored&signal_mask)/nnz(signal_mask));
    fprintf('percentage signal = %f \n', nnz(d_censored&signal_mask)/nnz(d_censored));
    fprintf('percentage crosstalk = %f \n', nnz(d_censored&ct_mask)/nnz(d_censored));
    fprintf('percentage background = %f \n', nnz(d_censored&noise_mask)/nnz(d_censored));

    [k_censored,t_censored] = get_kt(d_censored,im_size);

    % get ROM image only using detections left after background censoring
    % --> good starting point for depth estimation.
    rom_censored = get_rom(t_censored,3);
    % fill in ROM where data is missing (NaN pixels).
    rom_censored = roifill(rom_censored,imdilate(isnan(rom_censored),strel('disk',1,4)));
%     rom_censored(isnan(rom_censored)) = 0;

    % mean of detection times at each pixel --> input observations for
    % depth estimation. (times get added in NLL equation)
    t_censored_avg = nanmean(t_censored,3);
    t_censored_avg(isnan(t_censored_avg)) = 0;
    
    % adjust detection times for offset of pulse peak. 
    [~,pulse_peak] = max(pulse_hist);
    pulse_peak = pulse_bins(pulse_peak);
    t_censored_avg = t_censored_avg - 1 + (pulse_peak/time_res);
    t_censored_avg(t_censored_avg<0) = 0;
    
    %% 6. Initial depth estimation
    
    % max likelihood estimation
%     D_ML = t_censored_avg*time_res*c*0.5;
    
    % SPIRAL estimation
    D_AT = @(x) x; D_A = @(x) x;
    imag_out = SPIRALTAP(t_censored_avg,D_A,D1_tau, ...
        'noisetype',D_noisetype, ...
        'penalty',D_penalty, ...
        'maxiter',D1_maxiter,...
        'Initialization',rom_censored,... 
        'AT',D_AT,...
        'monotone',1,...
        'miniter',1,...
        'stopcriterion',3,...
        'tolerance',1e-8,...
        'alphainit',D1_ainit,...
        'alphaaccept',1e80,...
        'logepsilon',1e-10,...
        'saveobjective',1,...
        'savereconerror',1,...
        'savecputime',1,...
        'savesolutionpath',0,...
        'truth',D_truth,...
        'verbose',10);

    D_MAP = imag_out*time_res*c*0.5; % convert from time bins to depth
    figure;imshow(D_MAP,[95,102]);colorbar;colormap(jet)
    D_rmse = mean((D_MAP(:)-D_truth(:)).^2)^0.5
    
    %% Censor crosstalk (and background) detections.
    % based on likelihood that each detection is noise (due to crosstalk or background).

    sig_time = ceil(2*D_MAP./(c*time_res));
    is_sig_time = sparse(num_pixels,N);
    for n=1:N
        is_sig_time(:,n) = abs(detections(:,n)-sig_time(:))<=3;
    end

    signal_prob = I_MAP(:).*signal_level_im(:);
    noise_prob = noise_level_im(:)./max_timebin;
    lct = sparse(num_pixels,N);
    d_ct_censored = sparse(detections);
    for n=1:N
        lct(:,n) = (noise_prob + ct_prob(:,n))./(ct_prob(:,n)+(signal_prob .* is_sig_time(:,n)) + noise_prob);
        lct(isnan(lct(:,n)),n) = 0;
        d_ct_censored(:,n) = detections(:,n) .* (lct(:,n) < 0.5);
    end
    [k_ct_censored,t_ct_censored] = get_kt(d_ct_censored,im_size);
    ct_prob_im_ct_censored = reshape(full(sum(ct_prob.*(d_ct_censored~=0),2)),im_size);
    
    signal_ct_censored = nnz((~d_ct_censored)&signal_mask)/nnz(signal_mask);
    noise_ct_censored = nnz((~d_ct_censored)&noise_mask)/nnz(noise_mask);
    ct_ct_censored = nnz((~d_ct_censored)&ct_mask)/nnz(ct_mask);
    fprintf('accuracy = %f \n', nnz(((d_ct_censored~=0)==signal_mask)&detections)/nnz(detections));
    fprintf('fraction of background censored = %f \n', nnz((d_ct_censored==0)&noise_mask&d_censored)/nnz(noise_mask&d_censored));
    fprintf('fraction of crosstalk censored = %f \n', nnz((d_ct_censored==0)&ct_mask&d_censored)/nnz(ct_mask&d_censored));
    fprintf('fraction of signal retained = %f \n', nnz(d_ct_censored&signal_mask&d_censored)/nnz(signal_mask&d_censored));
    fprintf('percentage signal = %f \n', nnz(d_ct_censored&signal_mask)/nnz(d_ct_censored));
    fprintf('percentage crosstalk = %f \n', nnz(d_ct_censored&ct_mask)/nnz(d_ct_censored));
    fprintf('percentage background = %f \n', nnz(d_ct_censored&noise_mask)/nnz(d_ct_censored));
%     figure;hist(lct(signal_mask),100);
%     figure;hist(lct(ct_mask),100);
%     figure;hist(lct(noise_mask),100);

    %% Final intensity estimation

    I_ML = (log(N./(N-k_ct_censored))-ct_prob_im_ct_censored)./signal_level_im;
    
    I_MAP_ct_censored = SPIRALTAP_NEW(k_ct_censored,I_A,I2_tau,N, ...
        'noisetype',I_noisetype, ... % approx for binomial
        'penalty',I_penalty, ...
        'maxiter',I2_maxiter,... % `120
        'Initialization',I_ML,...
        'AT',I_AT,...
        'monotone',1,...
        'miniter',1,...
        'stopcriterion',3,...
        'tolerance',1e-10,...
        'alphainit',I2_ainit,...
        'alphaaccept',1e80,...
        'logepsilon',1e-10 + ct_prob_im_ct_censored,...
        'saveobjective',1,...
        'savereconerror',1,...
        'savecputime',1,...
        'savesolutionpath',0,...
        'truth',I_truth,...
        'verbose',10); 

    figure;imshow(I_MAP_ct_censored);

    I_ct_censored_psnr = psnr(I_MAP_ct_censored,I_truth)
    
    [counts,centers] = hist(I_MAP_ct_censored(:)-I_truth(:),100);
    [~,offset_ind] = max(counts);
    I_ct_censored_offset = centers(offset_ind);
    I_ct_censored_offset_psnr = psnr(I_MAP_ct_censored-I_ct_censored_offset,I_truth)
    
    %% Final depth estimation
    
    rom_ct_censored = get_rom(t_ct_censored,2);
%     rom_ct_censored(isnan(rom_ct_censored)) = 0;
    rom_ct_censored = roifill(rom_ct_censored,imdilate(isnan(rom_ct_censored),strel('disk',1,4)));
    
    t_ct_censored_avg = nanmean(t_ct_censored,3);
    t_ct_censored_avg(isnan(t_ct_censored_avg)) = 0;
    t_ct_censored_avg = t_ct_censored_avg - 1 + (pulse_peak/time_res);
    t_ct_censored_avg(t_ct_censored_avg<0) = 0;
    
    % SPIRAL estimation
    imag_out = SPIRALTAP(t_ct_censored_avg,D_A,D2_tau, ...
        'noisetype',D_noisetype, ...
        'penalty',D_penalty, ...
        'maxiter',D2_maxiter,...
        'Initialization',rom_ct_censored,... 
        'AT',D_AT,...
        'monotone',1,...
        'miniter',1,...
        'stopcriterion',3,...
        'tolerance',1e-8,...
        'alphainit',D2_ainit,...
        'alphaaccept',1e80,...
        'logepsilon',1e-10,...
        'saveobjective',1,...
        'savereconerror',1,...
        'savecputime',1,...
        'savesolutionpath',0,...
        'truth',D_truth,...
        'verbose',10);

    D_MAP_ct_censored = imag_out*time_res*c*0.5; % convert from time bins to depth
    figure;imshow(D_MAP_ct_censored,[95,102]);colorbar;colormap(jet)
    D_ct_censored_rmse = mean((D_MAP_ct_censored(:)-D_truth(:)).^2)^0.5
    
    save([filedir, sprintf('%d_results',iter)],'I_MAP','I_MAP_ct_censored','D_MAP','D_MAP_ct_censored','d_censored','d_ct_censored');
    save([filedir, sprintf('%d_performance',iter)],'I_psnr','I_ct_censored_psnr','D_rmse','D_ct_censored_rmse','I_offset','I_ct_censored_offset','I_offset_psnr','I_ct_censored_offset_psnr','signal_censored','noise_censored','ct_censored','signal_ct_censored','noise_ct_censored','ct_ct_censored');
end

