%combine highly correlated ICs to one
%used in twophoton_imageAnalysis

function [ data_corr, sm, sm_dsum] = finalMask_Jin(img_reg,mask_final, threshold)
%nmask = max(max(mask_final));
nmask = unique(mask_final);
nmask = nmask(nmask >0);
[npw, nph, ~] = size(img_reg);
%nmask = size(sel,2);
sm = zeros(npw,nph,length(nmask));
this_mask = 1;

%break up the 2D mask into a 3D mask
for ii=nmask %1:nmask
%     offsets = [1,-1, npw+1, npw-1,-npw+1,-npw-1, npw, -npw];
%     active_pixels = find(mask_final == ii);
%     active_pixels1 = bsxfun(@plus, active_pixels, offsets);
%     active_pixels1 = active_pixels1(:);
%     active_pixels1(mask_final(active_pixels1) < 0.5*ii)= [];
%     mask_final(active_pixels1) = ii;
    
    sm_placeholder = zeros(npw,nph);
    sm_placeholder(mask_final==ii) = 1;
    sm(:,:,this_mask)= sm_placeholder;
    this_mask = this_mask +1;
end

tc_avg = getTC(img_reg, sm, length(nmask));

data_corr = corrcoef(tc_avg);
%  figure; imagesc(data_corr);

%consolidate timecourses that are highly correlated
[i, j] = find(and(data_corr>threshold,data_corr<1)); %look for neurons which have sig corr but are not on the unity line
n = size(i,1);
sm_temp = reshape(sm,npw*nph,length(nmask));
remove_2 =[];
if n>1
    comp = [j(1); i(j == j(1))];
    total = unique(i); %its possible for a neuron to have sig corr with more than one other neuron. So find the total unique neurons with sig corrs.
    rest_ind = ismember(total, comp); %find all the masks not included in the comparison group
    rest_comp = total(~rest_ind); %
    
    %add all the masks in the comparison group to sm(:,:,comp(1))
    for k = 2:length(comp)
        sm_temp(find(sm(:,:,comp(k))),comp(1)) = 1; %put all of the masks from the comp group onto one mask
        sm(:, :, comp(1)) = reshape(sm_temp(:,comp(1)), npw ,nph); %update sm
    end
    
    if length(rest_comp) == 2
        sm_temp( find(sm(:,:,rest_comp(2))) ,rest_comp(1) ) = 1; %also add one of the rest comp masks to the all comp mask 
        sm(:, :, rest_comp(1)) = reshape(sm_temp(:,rest_comp(1)), npw ,nph); %update that sm mask
        remove_2 = rest_comp(2);
    elseif length(rest_comp) > 2
        for k = 1:length(rest_comp)
            if k > length(rest_comp)
                break
            else
                rest_comp1 = rest_comp(k);
                rest_comp2 = j(i == rest_comp1); %find all the masks which sig corr with rest_comp(k)
                rest_comp(ismember(rest_comp,rest_comp2))=[]; %removes rest_comp2 from rest_comp
                for kk = 1:length(rest_comp2)
                    sm_temp(find(sm(:,:,rest_comp2(kk))),rest_comp1) = 1; %put all of the masks from the comp group onto one mask
                    sm(:, :, rest_comp1) = reshape(sm_temp(:,rest_comp1), npw ,nph);
                end
                remove_2 = [remove_2; rest_comp2]; %keep track of 
            end
        end
    end
    remove_2 = [comp(2:end);remove_2]; 
    sm(:, :, remove_2) = [];
end
saveData = 0;
sm_dsum = plotMask(sm,saveData, 0,0);
end