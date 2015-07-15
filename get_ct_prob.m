function ct_prob = get_ct_prob(detections,p_ct,array_reps,im2array_inds)

ct_prob = sparse(size(detections,1),size(detections,2));
for n = 1:array_reps(1)
    for m = 1:array_reps(2)
        
        array_inds = im2array_inds(:, array_reps(1)*(m-1)+n);
        array_inds(array_inds==0) = [];
        
        ct_prob(array_inds,:) = compute_ct_prob(detections(array_inds,:), p_ct);
    end
end

end