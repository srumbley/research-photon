% Censor noise detections based on ROM and threshold. If the difference
% between a detection time at pixel (i,j) and rom(i,j) exceeds
% threshold(i,j), then the detection is censored as noise. 
% 
% Input t and output t_signal are size [height, width, k_max], where the
% first k(i,j) values at each pixel are the arrival times of detections,
% and the remaining values are NaN. 
% Output t_signal contains only detections that the threshold test
% determines to be signal. 

function d_censored = censor(detections,rom,threshold)

d_censored = sparse(detections);
for i=1:size(detections,2)
    inds = find(detections(:,i)~=0);
    
    d_censored(inds,i) = (abs(rom(inds)-nonzeros(detections(:,i))) < threshold(inds)) .* detections(inds,i);
end
% rom_diff = abs(bsxfun(@minus,rom(:),full(detections)) .* (detections ~= 0));
% is_signal = bsxfun(@lt,d_censored,threshold(:));
% d_censored = detections.*is_signal;

end