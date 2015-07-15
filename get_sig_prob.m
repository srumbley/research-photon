function sig_prob = get_sig_prob(detections, intensity, depth, signal_level_im, time_res)
sig_time = depth*2/(3e8*time_res);    
is_sig_time = sparse(detections);
for n=1:size(detections,2)
    is_sig_time(:,n) = abs(detections(:,n)-sig_time(:))<=2;
end
sig_prob = bsxfun(@times, intensity(:).*signal_level_im(:), is_sig_time);
end