function mask_final = processMask(mask_cell)
[npw, nph, nIC] = size(mask_cell);
mask_cell_temp = zeros(npw*nph, nIC);

% icasig = reshape(icasig,npw*nph,[]);
% ica_new = icasig;
% lr_offsets = [npw, -npw];
% ud_offsets = [1,-1, npw+1, npw-1,-npw+1,-npw-1];
%
% mask_cell_tt = zeros(npw*nph,nIC);
% icasig = reshape(icasig,npw*nph,[]);
% i = 1;
% while i <= nIC
% %     mask_temp   = zeros(npw*nph);
%
%     temp_ica = icasig(:,i);
%     peak_ica = max(temp_ica);
%     active_pixels = find(mask_cell(:,:,i)>0);
%
%     init = 1;
%     while ~isempty(init)
%
%         active_pixels1 = bsxfun(@plus, active_pixels, lr_offsets);
%
%         active_pixels1 = active_pixels1(:);
%
%         active_pixels1(temp_ica(active_pixels1) < 0.6*peak_ica) = [];
% %         active_pixels2 = bsxfun(@plus, active_pixels, ud_offsets);
% %         active_pixels2 = active_pixels2(:);
% %
% %         active_pixels2(temp_ica(active_pixels2) < 0.99*peak_ica) = [];
%
%         active_pixels1 = unique(active_pixels1);
%         active_pixels = [active_pixels(:); active_pixels1];
%         active_pixels(active_pixels <=0)=[];
%         active_pixels = unique(active_pixels);
%         init = [init;length(active_pixels1)];
%         if (length(init)>10)
%             if init(end) - init(end-1) < 4
%                 init=[];
%             end
%         end
%
%     end
%     if length(active_pixels) > 200
% %         mask_temp(active_pixels) = 1;
%         mask_cell_tt(active_pixels,i) = 1;
% %         mask_cell_tt(active_pixels2,i) = 1;
%         i = i + 1;
%     else
%         icasig(:,i) = [];
%         mask_cell(:,:,i) = [];
%         mask_cell_tt(:,i) = [];
%         nIC = nIC - 1;
%     end
%
% end
% mask_cell = reshape(mask_cell_tt,[npw,nph,nIC]);
for ic = 1:nIC
    if length(unique(reshape(mask_cell(:,:,ic),[1 npw*nph])))>2
        %         data_tc_temp = stackGetTimeCourses(img_reg,mask_cell(:,:,ic));
        %         data_corr_temp = corrcoef(data_tc_temp);
        ind_rem = 1:length(unique(reshape(mask_cell(:,:,ic),[1 npw*nph])))-1;
        for i = ind_rem
            ind_new = find(mask_cell(:,:,ic) == i);
            if length(ind_new) < 200
                continue
            else
                cat(3, mask_cell_temp, zeros(size(mask_cell_temp(:,:,1))));
                nIC = nIC+1;
                if i == 1
                    mask_cell_temp(find(mask_cell(:,:,ic)== i),ic) = 1;
%                     ica_new(:,ic) = icasig(:,ic);
                else
                    mask_cell_temp(find(mask_cell(:,:,ic)== ind_rem(i)),nIC) = 1;
%                     ica_new(:,nIC) = icasig(:,ic);
                end
            end
        end
    else
        mask_cell_temp(find(mask_cell(:,:,ic)),ic) = 1;
%         ica_new(:,ic) = icasig(:,ic);
    end
end
% lr_offsets = [npw, -npw];
% ud_offsets = [1,-1,npw-1, -npw+1, npw+1,-npw-1];
% 
% mask_cell_tt = zeros(npw*nph,nIC);
% 
% i = 1;
% while i <= nIC
% %     mask_temp   = zeros(npw*nph);
% 
%     temp_ica = icasig(:,i);
%     peak_ica = max(temp_ica);
%     active_pixels = find(mask_cell(:,:,i)>0);
% 
%     init = 1;
%     while ~isempty(init)
% 
%         active_pixels1 = bsxfun(@plus, active_pixels, lr_offsets);
% 
%         active_pixels1 = active_pixels1(:);
% 
%         active_pixels1(temp_ica(active_pixels1) < 0.55*peak_ica) = [];
%         active_pixels2 = bsxfun(@plus, active_pixels, ud_offsets);
%         active_pixels2 = active_pixels2(:);
% 
%         active_pixels2(temp_ica(active_pixels2) < 0.9*peak_ica) = [];
% 
%         active_pixels1 = unique(active_pixels1);
%         active_pixels = [active_pixels(:); active_pixels1;active_pixels2];
%         active_pixels(active_pixels <=0)=[];
%         active_pixels = unique(active_pixels);
%         init = [init;length(active_pixels1)];
%         if (length(init)>10)
%             if init(end) - init(end-1) < 4
%                 init=[];
%             end
%         end
% 
%     end
%     if length(active_pixels) > 200
% %         mask_temp(active_pixels) = 1;
%         mask_cell_tt(active_pixels,i) = 1;
% %         mask_cell_tt(active_pixels2,i) = 1;
%         i = i + 1;
%     else
%         icasig(:,i) = [];
%         mask_cell(:,:,i) = [];
%         mask_cell_tt(:,i) = [];
%         nIC = nIC - 1;
%     end
% 
% end
% mask_cell = reshape(mask_cell_tt,[npw,nph,nIC]);
%finds overlapping pixels of ROIs and based on correlations decides whether
%to group them or to split them- if splitting, then overlapping pixels are
%eliminated from both ROIs
nIC = size(mask_cell_temp,2);
mask_overlap = zeros(1,npw*nph);
mask_all = zeros(1,npw*nph);
min_size = 200;

for ic = 1:nIC
    
    ind_new = find(mask_cell_temp(:,ic))';
    if length(ind_new) >= min_size
        ind_old = find(mask_all);
        overlap = ismember(ind_old,ind_new);
        ind_both = find(overlap);
        if length(ind_both)>1
            ic_match = unique(mask_all(ind_old(ind_both)));
            for im = 1:length(ic_match)
                mask_all(ind_new) = ic;
                mask_all(ind_old(ind_both)) = 0;
                mask_overlap(ind_old(ind_both)) = 1;
            end
        else
            mask_all(ind_new) = ic;
        end
    end
end

% removes ICs smaller than 200 pixels, renumbers ROIs so in continuous ascending order
start = 1;
mask_final = zeros(size(mask_all));
for ic = 1:max(mask_all,[],2)
    ind = find(mask_all==ic);
    if length(ind)<200
        mask_overlap(mask_all==ic) = 1;
        mask_all(ind) = 0;
    end
    ind = find(mask_all==ic);
    if ~isempty(ind)
        mask_final(ind)=start;
        start= start+1;
    end
end
end