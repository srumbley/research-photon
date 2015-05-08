function [t_signal, t_signal_avg] = censor(t,rom,threshold)

rom_diff = abs(bsxfun(@minus,rom,t));
is_signal = bsxfun(@lt,rom_diff,threshold);
t_signal = t.*is_signal;
t_signal(t_signal==0) = NaN; % so that t_signal is NaN, not 0, for background detections

% has_data = sum(~isnan(t),3) >=1;
% has_data = repmat(has_data,1,1,size(t,3));
% t_signal(has_data==0) = NaN;
t_signal_avg = nanmean(t_signal,3);

% rom(isnan(rom)) = 0;
% num_pixels = size(t,1);
% t_signal = zeros(size(t));
% t_signal_avg = zeros(num_pixels);
% has_signal = zeros(num_pixels);
% % is_signal = false(size(t));
% for ii=1:num_pixels
%     if(mod(ii,100)==0)
%         fprintf([' row num = ' num2str(ii) '\n']);
%     end
%     
%     for jj=1:num_pixels
%         data_vec = t(ii,jj,:);
%         % 1. data missing from acquisition
%         if (sum(~isnan(data_vec))==0)
%             has_signal(ii,jj) = 0;
%         else
%             has_signal(ii,jj) = 1;
%         end
%         
%         % 2. censor data
%         if ( (rom(ii,jj) ~= 0) && (sum(~isnan(data_vec))==0) )
%             is_signal = abs(rom(ii,jj)-data_vec) < threshold(ii,jj);
%             signal_vec = data_vec.*is_signal;
%         else % rom is bad
%             has_signal(ii,jj) = 0;
%             signal_vec = zeros(1,1,length(data_vec));
%         end
%         
%         % 3. all data censored 
%         if (isempty(signal_vec))
%             has_signal(ii,jj) = 0;
%         end
%         t_signal(ii,jj,:) = signal_vec;
%         t_signal_avg(ii,jj) = mean(signal_vec);
%     end
% end
% 
% t_signal_avg(~has_signal) = NaN;
% sum(sum(sum((is_signal==(signal_mask.*(~isnan(t)))).*(~isnan(t)))))/sum(~isnan(t(:)))
end