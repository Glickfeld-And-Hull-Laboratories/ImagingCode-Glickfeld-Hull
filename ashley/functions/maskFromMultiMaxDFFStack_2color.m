function [mask_tagged, mask_other] = maskFromMultiMaxDFFStack_2color(maxDFF_stack)
%% interactive function to select cells from each stimulus driven max dF/F 
% in a y by x by n stimuli array. final output is the cell mask. 
% first image in maxDFF_stack is tagged cells (usually red channel)

xpix = size(maxDFF_stack,2);
ypix = size(maxDFF_stack,1);

tagImg = maxDFF_stack(:,:,1);
otherImg = maxDFF_stack(:,:,2:end);

nstim = size(otherImg,3);

% get tagged mask
bwout_tag = imCellEditInteractive(tagImg);
mask_tagged = bwlabel(bwout_tag);

% iterate through each max dF/F for each direction tuning stimulus
otherImg(repmat(bwout_tag,[1,1,nstim]) > 0) = 0;
% b = sum(imCellBuffer(tagImg,3),3);
% otherImg(b > 0) = 0;

bwout_other = zeros(ypix,xpix,nstim);
for istim = 1:nstim
    img_temp = otherImg(:,:,istim);
    if istim > 1
        prev_mask = bwout_other(:,:,istim-1); % set each previous cell == 0
        img_temp(prev_mask > 0) = 0;
        
        b = sum(imCellBuffer(prev_mask,3),3); % set 3 pixel buffer zone around each cell == 0
        img_temp(b > 0) = 0;
        
        bwout = imCellEditInteractive(img_temp); % select cells
        bwout_other(:,:,istim) = prev_mask+bwout;
        
    else
        bwout_other(:,:,istim) = imCellEditInteractive(img_temp); % select cells
    end
    
end
bwout_all = bwout_other(:,:,nstim);

if length(unique(bwout_all(:))) > 2
    error('you have overlap - set to zero')
end

mask_other = bwlabel(bwout_all);

end