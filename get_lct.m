function lct = get_lct(detections,params,intensity,depth)

if (rem(length(params),2)==1)
    error('Parameters should always go by pairs');
else
    for ii = 1:2:(length(params)-1)
        switch lower(params{ii})
            case 'p_ct';                p_ct                = params{ii+1}; 
            case 'signal_level_im';     signal_level_im     = params{ii+1}; 
            case 'noise_level_im';      noise_level_im      = params{ii+1};
            case 'time_res';            time_res            = params{ii+1};
            case 'pulse_timebins';      pulse_timebins      = params{ii+1}; 
            case 'max_timebin';         max_timebin         = params{ii+1};
            case 'im_size';             im_size             = params{ii+1};
            case 'array_size';          array_size          = params{ii+1};
            case 'array_reps';          array_reps          = params{ii+1};
            case 'im2array_inds';       im2array_inds       = params{ii+1};
                
        otherwise
                % Something wrong with the parameter string
                error(['Unrecognized option: ''', params{ii}, '''']);
        end
    end
end

ct_prob = get_ct_prob(detections, p_ct, array_reps, im2array_inds);
sig_prob = get_sig_prob(detections, intensity, depth, signal_level_im, time_res);
noise_prob = noise_level_im./max_timebin;
% pr_prob = bsxfun(@plus, sig_prob, noise_prob(:));
% lct = (noise_prob + ct_prob)./(ct_prob+pr_prob);
lct = sparse(detections);
for n=1:size(detections,2)
    lct(:,n) = ct_prob(:,n)./(ct_prob(:,n) + sig_prob(:,n) + noise_prob(:));
    lct(isnan(lct(:,n)),n) = 0;
end
end