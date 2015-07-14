filedir = 'results/full/tank2_2/';
results_files = dir([filedir, '*_results.mat']);
performance_files = dir([filedir, '*_performance.mat']);
n = length(results_files);

load([filedir, 'common_vars.mat']);
I_rmse_im = zeros(im_size);
D_rmse_im = zeros(im_size);
I_ct_censored_rmse_im = zeros(im_size);
D_ct_censored_rmse_im = zeros(im_size);
I_psnrs = zeros(n,1);
I_ct_censored_psnrs = zeros(n,1);
D_rmses = zeros(n,1);
D_ct_censored_rmses = zeros(n,1);

for i=1:n
    load([filedir, results_files(i).name]);
    I_rmse_im = I_rmse_im + (I_MAP - I_truth).^2;
    D_rmse_im = D_rmse_im + (D_MAP - D_truth).^2;
    I_ct_censored_rmse_im = I_ct_censored_rmse_im + (I_MAP_ct_censored - I_truth).^2;
    D_ct_censored_rmse_im = D_ct_censored_rmse_im + (D_MAP_ct_censored - D_truth).^2;
    load([filedir, performance_files(i).name]);
    I_psnrs(i) = I_psnr;
    I_ct_censored_psnrs(i) = I_ct_censored_psnr;
    D_rmses(i) = D_rmse;
    D_ct_censored_rmses(i) = D_ct_censored_rmse;
end

I_rmse_im = (I_rmse_im/n).^0.5;
D_rmse_im = (D_rmse_im/n).^0.5;
I_ct_censored_rmse_im = (I_ct_censored_rmse_im/n)./0.5;
D_ct_censored_rmse_im = (D_ct_censored_rmse_im/n)./0.5;

figure;imagesc(I_rmse_im);colorbar
figure;imagesc(D_rmse_im);colorbar
figure;imagesc(I_ct_censored_rmse_im);colorbar
figure;imagesc(D_ct_censored_rmse_im);colorbar
figure;plot(I_psnrs);hold on;plot(I_ct_censored_psnrs,'r')
figure;plot(D_rmses);hold on;plot(D_ct_censored_rmses,'r')
mean_I_psnr = mean(I_psnrs)
mean_D_rmse = mean(D_rmses)
mean_I_ct_censored_psnr = mean(I_ct_censored_psnrs)
mean_D_ct_censored_rmse = mean(D_ct_censored_rmses)

% function rmse_im = get_rmse_im(estimates,truth)
% 
% rmse = mean(bsxfun(@minus,estimates,truth).^2,3).^0.5;
% 
% end