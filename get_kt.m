function [k,t] = get_kt(detections, im_size)

k = reshape(full(sum(detections~=0,2)), im_size);
k_max = max(k(:));
t = sort(detections,2,'descend');
t = full(t(:,1:k_max));
t = reshape(t,im_size(1),im_size(2),k_max);
t(t==0) = NaN;

end