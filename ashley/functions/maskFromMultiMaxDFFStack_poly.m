function mask_cell = maskFromMultiMaxDFFStack_poly(maxDFF_stack)
%% interactive function to select cells from each stimulus driven max dF/F 
% in a y by x by n stimuli array. final output is the cell mask. 

xpix = size(maxDFF_stack,2);
ypix = size(maxDFF_stack,1);
nstim = size(maxDFF_stack,3);

% iterate through each max dF/F for each direction tuning stimulus
bwout_dir = zeros(ypix,xpix,nstim);

for istim = 1:nstim
    img_temp = maxDFF_stack(:,:,istim);
    if istim > 1
        prev_mask = bwout_dir(:,:,istim-1); % set each previous cell == 0
        img_temp(prev_mask > 0) = 0;
        
        b = sum(imCellBuffer(prev_mask,3),3); % set 3 pixel buffer zone around each cell == 0
        img_temp(b > 0) = 0;
        
        bwout = imCellPolyEditInteractive(img_temp); % select cells
        bwout_dir(:,:,istim) = prev_mask+bwout;
        
    else
        bwout_dir(:,:,istim) = imCellPolyEditInteractive(img_temp); % select cells
    end
    
end
bwout_all = bwout_dir(:,:,nstim);

if length(unique(bwout_all(:))) > 2
    error('you have overlap - set to zero')
end

mask_cell = bwlabel(bwout_all);

end