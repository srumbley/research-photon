% function rom = get_rom(t_cell,siz)
% tic;
% [nrows,ncols] = size(t_cell);
% 
% t_pad = [cell(nrows,siz), t_cell, cell(nrows,siz)];
% t_pad = [cell(siz,ncols+2*siz); t_pad; cell(siz,ncols+2*siz)];
% 
% rom = zeros(nrows,ncols);
% 
% for i=(1+siz):(nrows+siz)
%     if(mod(i,100)==0)
%         fprintf(['row=' num2str(i) '\n']);
%     end
%     
%     for j=(1+siz):(ncols+siz)
%         
%         t_neighborhood = t_pad((i-siz:i+siz), (j-siz:j+siz));
%         
%         raw_time = cell2mat(t_neighborhood(~cellfun('isempty',t_neighborhood(:))));
%         
%         filt_val = median(raw_time);
%         
%         if(isnan(filt_val))
%             rom(i-siz,j-siz) = 0;
%         else
%             rom(i-siz,j-siz) = filt_val;
%         end
%     end
% end
% toc;
% 
% end

% Get median of detection times within a sliding window of size [2*siz+1,2*siz+1].
% rom(i,j) is the median of the detection times in
% detections(i-siz:i+siz,j-siz:j+siz,:).

function rom = get_rom(t,siz)
tic;
[nrows, ncols, k_max] = size(t);

d_pad = nan(nrows+2*siz, ncols+2*siz, k_max);
d_pad(siz+1:end-siz, siz+1:end-siz, :) = t;

rom = zeros(nrows,ncols);

for i=(1+siz):(nrows+siz)    
    if(mod(i,100)==0)
        fprintf(['row = ' num2str(i) '\n'])
    end

    for j=(1+siz):(ncols+siz)
        if (j==4)
            if (i==278)
            
                disp('asdf');
            end
        end
        % get all detections in neighborhood
        d_neighborhood = d_pad((i-siz:i+siz),(j-siz:j+siz),:);
        
        % get median of all detection times in neighborhood
        rom(i-siz,j-siz) = nanmedian(d_neighborhood(:));
        
    end
end
% 
% toc;
% 
% % using pixel index in order to use parfor. only slightly faster...
% % parfor p=1:(nrows*ncols)
% %     % get image subscripts - this is the left corner of the neighborhood
% %     % window in the padded image
% %     [i,j] = ind2sub([nrows,ncols],p);
% %     
% %     % get detections in neighborhood
% %     d_neighborhood = d_pad((i:i+2*siz),(j:j+2*siz),:);
% %     
% %     % get median of all detections in neighborhood
% %     rom(p) = nanmedian(d_neighborhood(:));
% % end
% end
% % 
% % this is slower :( but so much prettier...
% % function rom = get_rom(t, siz)
% % tic;
% % rom = blockproc(t,[1,1],@(block_struct) nanmedian(block_struct.data(:)), ...
% %                 'BorderSize',[siz,siz],...
% %                 'TrimBorder',false,...
% %                 'PadMethod',NaN,...
% %                 'UseParallel',true);
% % toc;
% % end
