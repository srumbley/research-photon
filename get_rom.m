% Get median of detection times within a sliding window of size [2*siz+1,2*siz+1].
% rom(i,j) is the median of the detection times in
% detections(i-siz:i+siz,j-siz:j+siz,:).

function rom = get_rom(t,siz)

tic;
[nrows, ncols, k_max] = size(t);

% ...tried using weighted median - results not that good.
% weights = [1:siz+1,siz:-1:1];
% weights = weights'*weights;
% weights = repmat(weights,1,1,k_max);

d_pad = nan(nrows+2*siz, ncols+2*siz, k_max);
d_pad(siz+1:end-siz, siz+1:end-siz, :) = t;

rom = zeros(nrows,ncols);
% local_var = zeros(nrows,ncols);

for i=(1+siz):(nrows+siz)    
    if(mod(i,100)==0)
        fprintf(['row = ' num2str(i) '\n'])
    end

    for j=(1+siz):(ncols+siz)
        % get all detections in neighborhood
        d_neighborhood = d_pad((i-siz:i+siz),(j-siz:j+siz),:);
        d_neighborhood(siz+1,siz+1,:) = NaN(1,1,k_max);
%         inds = ~isnan(d_neighborhood);
%         if sum(inds(:))>0
%             rom(i-siz,j-siz) = weighted_median(d_neighborhood(inds), weights(inds));
%         end
        % get median of all detection times in neighborhood
        rom(i-siz,j-siz) = median(d_neighborhood(~isnan(d_neighborhood)));
%         local_var(i-siz,j-siz) = nanvar(d_neighborhood(:));
    end
end

toc;

% using pixel index in order to use parfor. only slightly faster...
% parfor p=1:(nrows*ncols)
%     % get image subscripts - this is the left corner of the neighborhood
%     % window in the padded image
%     [i,j] = ind2sub([nrows,ncols],p);
%     
%     % get detections in neighborhood
%     d_neighborhood = d_pad((i:i+2*siz),(j:j+2*siz),:);
%     
%     % get median of all detections in neighborhood
%     rom(p) = nanmedian(d_neighborhood(:));
% end
end
