% Censor noise detections based on ROM and threshold. If the difference
% between a detection time at pixel (i,j) and rom(i,j) exceeds
% threshold(i,j), then the detection is censored as noise. 
% 
% Input t and output t_signal are size [height, width, k_max], where the
% first k(i,j) values at each pixel are the arrival times of detections,
% and the remaining values are NaN. 
% Output t_signal contains only detections that the threshold test
% determines to be signal. 

function t_signal = censor(t,rom,threshold)

rom_diff = abs(bsxfun(@minus,rom,t));
is_signal = bsxfun(@lt,rom_diff,threshold);
t_signal = t.*is_signal;
t_signal(t_signal==0) = NaN; % so that t_signal is NaN, not 0, for censored detections

end