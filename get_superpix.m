% Group detections into superpixels.
%
% Inputs:
% - detections: output of simulate_detections. sparse matrix of detection
% times (size = num_pixels x num_pulses).
% - im_size: size of original image (height x width).
% - superpix_size: [superpix_height, superpix_width].
%
% Outputs: analogous to outputs of simulate_detections, but at reduced
% image size (im_size(1)/superpix_size(1), im_size(2)/superpix_size(2)). 

function [k,t,superpix_detections] = get_superpix(detections,im_size,superpix_size)
tic;
num_pulses = size(detections,2);
size_out = ceil(im_size./superpix_size);
superpix_detections = spalloc(size_out(1)*size_out(2), num_pulses, nnz(detections));

% pre-compute superpixel indices
inds = reshape((1:im_size(1)*im_size(2)),im_size);
block_inds = im2col(inds,superpix_size,'distinct');

% subtract max+1 from nonzero values --> find min, add max+1 back
max_t = max(detections(:)) + 1;
detections2 = spfun(@(x) x-max_t, detections);

for p = 1:num_pulses
    if(mod(p,100)==0)
        fprintf('superpixing pulse #%d \n',p);
    end
    
    detections_p = detections2(:,p);
    detections_cols = detections_p(block_inds);
    superpix_detections(:,p) = min(detections_cols)';
end
superpix_detections = spfun(@(x) x + max_t, superpix_detections);
toc;

k = full(sum(superpix_detections > 0, 2));
k_max = max(k);
k = reshape(k, size_out);

t = sort(superpix_detections,2,'descend');
t = full(t(:,1:k_max));
t(t==0) = NaN;
t = reshape(t, size_out(1), size_out(2), k_max);
end