%% 1. load ground truth data
load('ground_truth_man_intensity.mat'); % I_truth
load('gt_depth_man.mat'); % D_truth
I_truth = I_truth(1:2:end,1:2:end);
D_truth = D_truth(1:2:end,1:2:end);
I_truth = I_truth-min(I_truth(:));
I_truth = I_truth/max(I_truth(:));

im_size = size(I_truth);

%% 2. set simulation parameters

% if array size is smaller than image size
array_size = [64,64];
array_reps = ceil(im_size./array_size);

N = 1000;
% q = 0.17;
q_array = 0.17 + randn(array_size)*0.05^2;
q = repmat(q_array,array_reps);
q = q(1:im_size(1),1:im_size(2));

h_array = fspecial('gaussian',array_size,0.75*array_size(1)); % Gaussian beam falloff
h_array = h_array/max(h_array(:));
h = repmat(h_array,array_reps);
h = h(1:im_size(1),1:im_size(2));
S = 0.09 * h;

b = 0.001;
d = 0;

signal_level = q.*S;
noise_level = q.*b + d;

Tr = 100e-9;
time_res = 8e-12;
pulse_args = {'pulse_gauss', 270e-12, [-500e-12,500e-12]};

deadpix = rand(array_size);
deadpix = deadpix >= 0.001;
deadpix = repmat(deadpix,array_reps);
deadpix = deadpix(1:im_size(1),1:im_size(2));

%% 3. generate simulated data
[k,t,detections,signal_mask] = simulate_detections(I_truth,D_truth,N,signal_level,noise_level,Tr,time_res,pulse_args);

% kill dead pixels
k = k.*deadpix;
t = bsxfun(@times,t,deadpix);
signal_mask = bsxfun(@times,signal_mask,deadpix(:));

fprintf(['fraction of pixels with no detections: ' num2str(sum(sum(k==0))/length(k(:))) '\n']);
fprintf(['max # detections: ' num2str(max(k(:))) '\n']);
fprintf(['mean # detections: ' num2str(mean(k(:))) '\n']);

%% 4. estimate intensity

fprintf('\n Estimating intensity... \n');

% max likelihood estimation
I_ML = (k./(N-k)-noise_level)./(q.*S);

% SPIRAL parameters
tau = 7; ainit = 100; max_iter = 30; 
pen_type = 'tv'; noisetype = 'poisson';

% if using wavelets for smoothing, input these params to SPIRALTAP
% wav = daubcqf(2);
% W = @(x) midwt(x,wav);
% WT = @(x) mdwt(x,wav);

A = @(x) signal_level.*x + noise_level; % observation matrix
AT = @(x) signal_level.*x; 

tic;

I_MAP = SPIRALTAP_NEW(k,A,tau,N, ...
    'noisetype',noisetype, ... % approx for binomial
    'penalty',pen_type, ...
    'maxiter',max_iter,... % `120
    'Initialization',I_ML./max(I_ML(:)),...
    'AT',AT,...
    'monotone',1,...
    'miniter',1,...
    'stopcriterion',3,...
    'tolerance',1e-10,...
    'alphainit',ainit,...
    'alphaaccept',1e80,...
    'logepsilon',1e-10,...
    'saveobjective',1,...
    'savereconerror',1,...
    'savecputime',1,...
    'savesolutionpath',0,...
    'truth',zeros(im_size),...
    'verbose',5,...     % print info every x iterations
    'showfigiter',5);   % show figure every x iterations
%     'savetlv',1
toc;

figure;imshow(I_MAP,[]);colorbar

%% 5. ROM filtering

I_MAP_scaled = I_MAP/max(I_MAP(:));
threshold = 50./(I_MAP_scaled+0.5);
rom = get_rom(t,3);
t_signal = censor(t,rom,threshold);
t_signal_avg = nanmean(t_signal,3);

%% 6. Depth estimation

% max likelihood estimation
D_ML = t_signal_avg*8e-12*3e8*0.5;

% SPIRAL estimation
rom(isnan(rom)) = 0;
t_signal_avg(isnan(t_signal_avg)) = 0;

tau = 20; ainit = 0.05; max_iter = 50;

set_penalty = 'tv';
AT = @(x) x; A = @(x) x;
tic;
imag_out = SPIRALTAP(t_signal_avg,A,tau, ...
    'noisetype','gaussian', ...
    'penalty',set_penalty, ...
    'maxiter',max_iter,...
    'Initialization',rom,... 
    'AT',AT,...
    'monotone',1,...
    'miniter',1,...
    'stopcriterion',3,...
    'tolerance',1e-12,...
    'alphainit',ainit,...
    'alphaaccept',1e80,...
    'logepsilon',1e-10,...
    'saveobjective',1,...
    'savereconerror',1,...
    'savecputime',1,...
    'savesolutionpath',0,...
    'truth',zeros(im_size),...
    'verbose',5);
toc;
D_MAP = imag_out*8e-12*3e8*0.5;
% figure;imshow(D_MAP,[3.68,3.82]);colorbar % mug
figure;imshow(D_MAP,[4.25,4.6]);colorbar % man
