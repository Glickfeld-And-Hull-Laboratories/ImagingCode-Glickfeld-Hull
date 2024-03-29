% main: process movie as following:
% 1. registeration -- find the most stable frame out of 100
% 2. PCA -- #PCA depends on how dense dendrite population, 100-200
% 3. ICA -- #IC is less than #PCA and with 30HZ, the temperal-spatial ratio
% is set to be ~1
% 4. Apply gaussian filter to ICA
% 5. Threshold ICA signal to get binary mask
% 6. Process mask (break multiple components on single layer to multiple,
% combine highly correlated ones and remove small ones)
% 7. Extrack timecourse (manually inspect TC and remove bad ones) and create a singal layer mask for display

clear
file_info_WF;
img_dir = 'Z:\home\jake\Data\WidefieldImaging\GCaMP\';

for sub = [1]%size(mouseID,2)
    for rID = 1
        file_info_WF;
        usFacs = 100;
        out_dir  = fullfile('C:','Users','ziye','Documents','MATLAB','WF_Analysis',[date{sub} '_', mouseID{sub}],'\');
        
        mouse = mouseID{sub};
        dateID = date;
        date = dateID{sub};
        img_dir = dir([img_dir date '_' mouse]);
        img = readtiff(img_dir, 'uint16');
%         [img, skip_run, img_fn] = loadFile(sub, rID);
        
        if skip_run == 1
            continue
        else
            if ~exist(out_dir)
                mkdir(out_dir);
            end
            
            if exist([out_dir, 'ROI_xy.mat'],'file') == 2
                load([out_dir, 'ROI_xy.mat']);
            else
                [ROI_x, ROI_y] = get_2P_ROI(img); % get the ROI -  must be a rectangle
                save([out_dir 'ROI_xy.mat'],  'ROI_x', 'ROI_y');
            end
            %
            img = img(ROI_x,ROI_y,:);
            %             img = img(:,:,22976:end);
            %%
            
            [img_reg, pre] = homomorphicFilter(img);
            [img_reg] = spatialFilter(img_reg);
            
            if exist([out_dir, 'Reg_out.mat'],'file') == 2
                load([out_dir, 'Reg_out.mat']);
                [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,reg_out);
            else
                img = double(img);
                rf = 30*100;
                ref30 = img(:,:,10:100:rf);
                sf = 100*100+10;
                samp100 = img(:,:,11:100:sf);
                dshift = [];
                for r = 1:size(ref30,3)
                    
                    [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
                    dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
                    
                end
                
                %%
                
                min_f = find(dshift == min(dshift));
                img_ref = ref30(:,:,min_f);
                [reg_out, img_reg] = stackRegister(img, img_ref);
                save([out_dir, 'Reg_out.mat'], 'reg_out','img_ref');
            end
            
            clear img
            reg_sum = sum(img_reg,3);
%             dis = sqrt((reg_out(:,3).^2)+(reg_out(:,4).^2));
%             img_red = img_reg(:,:,dis<1.5);
            
            [npw, nph, nt] = size(img_reg);
            
            nPCA = 100; %100
            img_pca = img_reg(:,:,1:1:end);
            nf = size(img_pca,3);
            [mixedsig, mixedfilters, CovEvals, ~, ~, ...
                ~] = CellsortPCA2(img_pca,[1 nf], nPCA,[], out_dir, mouseID{sub}, []);
%                          [mixedsig, mixedfilters, CovEvals, ~, ~, ~] = CellsortPCA2(img_reg,[1 nt], nPCA, [2 1], out_dir, mouseID{sub}, []);
            
            PCuse = 1:nPCA;
            %             [PCuse] = CellsortChoosePCs(mixedfilters);
            
            %% Compute independent components
            mm = 0.8;
            nIC = 50;
            termtol = 0.00001;
            maxrounds = 1000;
            
            [ica_sig, mixedfilters, ica_A, numiter] = CellsortICA(mixedsig, ...
                mixedfilters, CovEvals, PCuse, mm, nIC, [], termtol, maxrounds);
            
            icasig = permute(mixedfilters,[2,3,1]);
            
%             save([out_dir, 'ICA.mat'], 'ica_sig', 'icasig'); %clear ica_sig;
            
            nIC = size(icasig,3);
            
            icasig = stackFilter(icasig);
            
            mask_cell = zeros(size(icasig));
            sm_logical = zeros(npw,nph);
            cluster_threshold = 90; %94.5
            
            for ic = 1:nIC
                icasig(:,:,ic) = imclearborder(icasig(:,:,ic));
                sm_logical((icasig(:,:,ic)>mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=1;
                sm_logical((icasig(:,:,ic)<=mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=0;
                sm_logical = logical(sm_logical);
                mask_cell(:,:,ic) = bwlabel(sm_logical);
            end
            
            %                         mask_cell = zeros(size(icasig));
            %                         for ic = 1:nIC
            %                             bwimgcell = imCellEditInteractive(icasig(:,:,ic),[]);
            %                             mask_cell(:,:,ic) = bwlabel(bwimgcell);
            %                             close all
            %                         end
            mask_final = processMask(mask_cell);
            
            mask_raw = reshape(mask_final, npw, nph);
            figure;imagesc(mask_raw);truesize
            %             h = impoly;
            %             junk_mask = (createMask(h));
            %             mask_raw(junk_mask) = 0;
            %             figure;imshow(mat2gray(kk));truesize
            %             hold all; h = contour(mask_raw);
            %             set(h, 'AlphaData',0.3);
            %
            %             saveas(h, [out_dir, 'mask_overlay.fig']);
            %             print([out_dir, 'mask_overlay.eps'],'-depsc')
            %             figure; h = imagesc(mask_raw);truesize
            %             saveas(h, [out_dir, 'mask.fig']);
            %             print([out_dir, 'mask.eps'],'-depsc')
            %
            %             mask_final = reshape(mask_raw, 1, npw*nph);
            %             mask_raw(mask_raw == 34)
            %            mask_final = reshape(mask_raw, 1, npw*nph);
            threshold = 0.8;
            
            [ ~, mask3D, ~] = finalMask(img_reg(:,:,1:10:end), mask_final, threshold, out_dir);
            
            %% Plotting TCs
            nmask = size(mask3D,3);
            FrameRate = 30;
            tc_avg = getTC(img_reg, mask3D, nmask);
            
            %             ICbad = [];
            %             saveData = 0;
            %             kl = floormask/5);
            %             for k = 1:kl
            %                 plotTC(tc_avg, mask3D, reg_sum, (k-1)*5+1:k*5, FrameRate, out_dir,saveData);
            %                 ICbad_input = input('Number of bad IC ', 's');
            %                 ICbad = [ICbad str2num(ICbad_input)];
            %
            %             end
            %             if mod(nmask,5) ~= 0
            %                plotTC(tc_avg, mask3D, reg_sum, 5*kl+1:nmask, FrameRate, out_dir,saveData);
            %                ICbad_input = input('Number of bad IC ', 's');
            %                ICbad = [ICbad str2num(ICbad_input)];
            %             end
            %             close all
            %             tc_avg(:,ICbad) = []; mask3D(:,:,ICbad) = [];
            saveData = 1;
            plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
            mask_flat = plotMask(mask3D, saveData, out_dir);
            
            data_corr = corrcoef(tc_avg);
            figure; fig = imagesc(data_corr);
            saveas(fig, [out_dir, 'data_corr.fig']);
            print([out_dir, 'data_corr.eps'],'-depsc')
            
            mask_final = processMask(mask3D);
            sz = [npw, nph];
            save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
            save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
            close all
            
            %             sub = 16;
            %             rID = 1;
            mouse = mouseID{sub};
            dateID = date;
            date = dateID{sub};
            data_dir = fullfile('\\crash\data\home\jake\Analysis\2P Analysis\Lever');
            subfold = fullfile(data_dir,[dateID{sub} '_' mouseID{sub} '_run' runID{rID} '\']);
            % foldName = subfold.name;
            % data_dir = fullfile(data_dir,foldName,'\');
            dataName   = dir(fullfile(subfold,'*frame_times.mat'));
            load([subfold, dataName.name]);
            
            tc_dir  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
            dest = tc_dir;
            dest_sub = dest;
            %             load([tc_dir, 'ROI_TCs.mat']);
            
            getTC_events;
            HAD_2P_TC_quantification;
            close all
            clear
        end
    end
end

% %             smwidth = 3;
% %             thresh = 6;
% %             arealims = 200;
% %             plotting = 1;
% %             weight = 0;
% %
% %             [ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(mixedfilters, smwidth, thresh, arealims, plotting, weight);
% % %
% %             subtractmean = 0; normal = 0;
% %
% %             cell_sig = CellsortApplyFilter(img_reg(:,:,1:4:end), ica_segments, [], movm, subtractmean, normal);
% %
% %             [data_corr, sm, mask_flat] = finalMask2( cell_sig', ica_segments, threshold, out_dir);
% %
% % %             ica_seg = permute(sm, [3,1,2]);
% % %
% % %             cell_sig = CellsortApplyFilter(img_reg, ica_seg, [], movm, subtractmean);
% %
% %             deconvtau = 0;
% %             spike_thresh = 2;
% %             normalization = 0;
% %             dt = 0.1;
% %             [~, spt, spc] = CellsortFindspikes(cell_sig, spike_thresh, dt, deconvtau, normalization);
% %
% %             %% Show results
% %             ICbad = [];
% %             for i = 1:nIC/5
% %                 figure;
% %                 CellsortICAplot('contour', ica_seg, cell_sig, movm, [], dt, 1, 2, [(i-1)*5+1:i*5], spt, spc);
% %                 ICbad_input = input('Number of bad IC ', 's');
% %                 ICbad = [ICbad str2num(ICbad_input)];
% %             end
% %             close all
% %
% %             ica_seg(ICbad2,:,:) = []; cell_sig(ICbad2,:) = [];
% %             [~, spt, spc] = CellsortFindspikes(cell_sig, spike_thresh, dt, deconvtau, normalization);
% %             figure;CellsortICAplot('contour', ica_seg, cell_sig, movm, [], dt, 1, 1, 1:size(ica_seg,1) , spt, spc);
% %
% % %             save([out_dir, 'Results_new.mat'], 'cell_sig', 'ica_seg', 'mask_flat', 'data_corr');
% %
% %
% %
% %             cell_sig = cell_sig';
