%% get path names
close all;clear all;clc;

ds = 'AdaptPhaseRev_ExptList';
eval(ds)
nexp = length(expt);
rc = behavConstsAV;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

iexp = 5;
             %%
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        area = expt(iexp).img_loc{1};
        ImgFolder = expt(iexp).prFolder;
        adaptFolder = expt(iexp).adaptFolder;
        time = expt(iexp).prTime;
        nrun = length(ImgFolder);
        frameRateHz = params.frameRate;

        run_str = catRunName(ImgFolder, nrun);
        ad_run_str = catRunName(adaptFolder, size(adaptFolder,2));

            fprintf(['2P imaging Phase analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
            for irun=1:nrun
                fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
            end

                %% load
                tic
                data = [];
                clear temp
                trial_n = [];
                offset = 0;
                for irun = 1:nrun
                    %CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder{irun}];
                    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
                    cd(CD);
                    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
                    load(imgMatFile);
                    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
                    load(fName);
                    
                    temp(irun) = input;
                    nOn = temp(irun).nScansOn;
                    nOff = temp(irun).nScansOff;
                    ntrials = size(temp(irun).tGratingDirectionDeg,2);
                    nframes = ntrials*(nOn+nOff);


                    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
                    data_temp = sbxread(imgMatFile(1,1:11),0,nframes);
                    if size(data_temp,1)== 2
                        data_temp = data_temp(1,:,:,:);
                    end
                    data_temp = squeeze(data_temp);
                    data = cat(3,data,data_temp);
                    fprintf('Complete')
                end
                input = concatenateDataBlocks(temp);
                fprintf('\nAll runs read\n')
                fprintf([num2str(size(data,3)) ' total frames\n'])
                clear data_temp
                clear temp

                toc

                %% register to adapt experiment

                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ad_run_str], [date '_' mouse '_' ad_run_str '_reg_shifts.mat']))
                [out, data_reg] = stackRegister(data,data_avg);
                data_reg_avg = mean(data_reg,3);
                mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
                save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'data_reg_avg')
                save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
                clear data
                %% test stability
                
                figure; 
                movegui('center')
                subplot(2,2,1);
                imagesc(data_reg_avg);
                title('PhaseRev run avg')
                subplot(2,2,2);
                imagesc(data_avg)
                title('Adapt run avg')
                sz = size(data_avg);
                rgb = zeros(sz(1),sz(2),3);
                rgb(:,:,1) = data_reg_avg./max(data_reg_avg(:));
                rgb(:,:,2) = data_avg./max(data_avg(:));
                subplot(2,2,3);
                image(rgb)
                title('Overlay')

                print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

                %% use adapt mask to get TCs

                fprintf(['Loading masks from adapt runs: ' cell2mat(adaptFolder) '\n'])

                % loads 'mask_cell', 'mask_np'
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ad_run_str], [date '_' mouse '_' ad_run_str '_mask_cell.mat']))
                fprintf('Cell and neuropil masks loaded\n')

                nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
                fprintf([num2str(nCells) ' total cells selected\n'])
                fprintf('Cell segmentation complete\n')

                % neuropil subtraction
                down = 5;
                sz = size(data_reg);

                data_tc = stackGetTimeCourses(data_reg, mask_cell);
                data_reg_down  = stackGroupProject(data_reg,down);
                data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
                nCells = size(data_tc,2);
                np_tc = zeros(sz(3),nCells);
                np_tc_down = zeros(floor(sz(3)./down), nCells);
                for i = 1:nCells
                     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
                     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
                     fprintf(['Cell #' num2str(i) '%s/n']) 
                end
                %get weights by maximizing skew
                ii= 0.01:0.01:1;
                x = zeros(length(ii), nCells);
                for i = 1:100
                    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
                end
                [max_skew ind] =  max(x,[],1);
                np_w = 0.01*ind;
                npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
                clear data_reg data_reg_down

                save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

                fprintf('\nNeuropil subtraction complete\n')

                clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

            %% Phase reversal analysis
            nOn = double(input.nScansOn);
            nOff = double(input.nScansOff);
            phaseCyc = double(input.nScansPhaseCyc);
            ntrials = length(input.tGratingDirectionDeg);
            nCells = size(npSub_tc,2);

            cycPerTrial = floor(nOn/(phaseCyc*2));
            data_trial = permute(reshape(npSub_tc,[nOn+nOff ntrials nCells]),[1 3 2]);
            data_f = mean(data_trial(nOff./2:nOff,:,:),1);
            data_dfof = (data_trial-data_f)./data_f;

            data_dfof_cyc = zeros((phaseCyc.*2)+phaseCyc/2, nCells, ntrials, cycPerTrial);
            for icyc = 1:cycPerTrial
                data_dfof_cyc(:,:,:,icyc) = data_dfof(nOff+(phaseCyc/2)+((icyc-1).*(phaseCyc.*2)):nOff+phaseCyc+(icyc.*(phaseCyc.*2))-1,:,:);
            end
            data_dfof_cycavg = mean(data_dfof_cyc,4);

            dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
            dirs = unique(dir_mat);
            nDir = length(dirs);

            tStimNum = celleqel2mat_padded(input.tStimulusNumber);
            nPhase = input.gratingStartingPhaseStepN;
            phase_mat = zeros(1,ntrials);
            for itrial = 1:ntrials
                if tStimNum(itrial) < (nDir.*nPhase)
                    temp = tStimNum(itrial);
                else
                    temp = mod(tStimNum(itrial),nDir.*nPhase);
                end
                if temp < nPhase
                    phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*temp);
                else
                    phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*(mod(temp,input.gratingStartingPhaseStepN)));
                end
            end
            phases = unique(phase_mat);

            base_win = cell(1,2);
            resp_win = cell(1,2);
            base_win{1} = [phaseCyc/2-2:phaseCyc/2+3];
            resp_win{1} = [(phaseCyc/2)+ceil(frameRateHz/3):(phaseCyc/2)+ceil(frameRateHz/2)];
            base_win{2} = [1.5*phaseCyc-2:1.5*phaseCyc+3];
            resp_win{2} = [1.5*phaseCyc+ceil(frameRateHz/3):1.5*phaseCyc+ceil(frameRateHz/2)];
            tt = (1-phaseCyc/2:phaseCyc*2).*(1000/frameRateHz);

            data_dfof_phasedir = zeros(2.5.*phaseCyc, nCells, nPhase, nDir);
            h_dir = zeros(nCells,nPhase,nDir,2);
            p_dir = nan(nCells,nPhase,nDir,2);
            phasedir_resp_avg = zeros(nCells,nPhase,nDir,2,2);
            trial_n = zeros(nPhase,nDir);
            for iPhase = 1:nPhase
                ind_phase = find(phase_mat == phases(iPhase));
                for iDir = 1:nDir
                    ind_dir = find(dir_mat == dirs(iDir));
                    ind = intersect(ind_phase,ind_dir);
                    trial_n(iPhase,iDir) = length(ind);
                    data_dfof_phasedir(:,:,iPhase,iDir) = mean(data_dfof_cycavg(:,:,ind)-mean(data_dfof_cycavg(base_win{1},:,ind),1),3);
                    for i = 1:2
                        phasedir_resp_avg(:,iPhase,iDir,i,1) = squeeze(mean((mean(data_dfof_cycavg(resp_win{i},:,ind),1)-mean(data_dfof_cycavg(base_win{i},:,ind),1)),3));
                        phasedir_resp_avg(:,iPhase,iDir,i,2) = squeeze(std((mean(data_dfof_cycavg(resp_win{i},:,ind),1)-mean(data_dfof_cycavg(base_win{i},:,ind),1)),[],3))./sqrt(length(ind));
                        if length(ind)>3
                            if nDir>1
                                alpha = 0.05./nDir-1;
                            else
                                alpha = 0.05;
                            end
                            [h_dir(:,iPhase,iDir,i), p_dir(:,iPhase,iDir,i)] = ttest2(squeeze(mean(data_dfof_cycavg(resp_win{i},:,ind),1)), squeeze(mean(data_dfof_cycavg(base_win{i},:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', alpha);
                        end
                    end
                end
            end

            [max_val, max_ind] = max(max(max(phasedir_resp_avg(:,:,:,1),[],4),[],2),[],3);

            resp_ind_phase = find(sum(sum(sum(h_dir,2),3),4));
            [n n2] = subplotn(length(resp_ind_phase));
            figure;
            for iCell = 1:length(resp_ind_phase)
                subplot(n,n2,iCell)
                temp = reshape(permute(squeeze(h_dir(resp_ind_phase(iCell),:,:,:)),[1 3 2]), [nPhase nDir.*2]);
                imagesc(temp)
                set(gca, 'YTick', 1:nPhase, 'YTickLabel', num2str(phases'))
                ylabel('Phase')
                set(gca, 'XTick', 1.5:2:7.5, 'XTickLabel', num2str(dirs'))
                title([num2str(resp_ind_phase(iCell))])
            end
            suptitle([date ' ' mouse '- Phase reversal'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sigRespByPhase&Dir.pdf']),'-dpdf','-bestfit');

            save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof', 'resp_win', 'base_win', 'tt', 'h_dir', 'p_dir', 'resp_ind_phase', 'data_dfof_phasedir', 'phasedir_resp_avg', 'max_ind', 'nCells')
            save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'phase_mat', 'phases', 'nPhase', 'dir_mat', 'dirs', 'nDir', 'nOn', 'nOff','frameRateHz', 'phaseCyc')

            %% measure F2/F1
            f1 = zeros(nDir,nCells);
            f2 = zeros(nDir,nCells);
            for iCell = 1:nCells
                for iDir = 1:nDir
                    cyc = squeeze(data_dfof_phasedir(phaseCyc./2+1:end,iCell,:,iDir))';
                    [f1(iDir,iCell),f2(iDir,iCell),f1ang,projectedf1,f1mat,f2mat] = compcontrastrevf1f2(cyc);
                end
            end
            [max_val_f1 max_dir_f1] = max(f1',[],2);
            
            figure;
            movegui('center')
            [n n2] = subplotn(nCells);
            for iCell = 1:nCells; 
                subplot(n, n2, iCell)
                plot(tt, squeeze(data_dfof_phasedir(:,iCell,:,max_dir_f1(iCell))))
                title(num2str(chop(max_val_f1(iCell)./std(data_dfof_phasedir(1:20,iCell,1,1),[],1),2)))
                if find(resp_ind_phase==iCell)
                    vline(0)
                end
            end
            suptitle([date ' ' mouse '- Phase reversal'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgRespByPhase_prefDir.pdf']),'-dpdf','-fillpage');

            
            
            ind_f1 = sub2ind(size(f1'),(1:size(f1',1))',max_dir_f1);
            f1_inv = f1';
            f2_inv = f2';
            f1_dir = f1_inv(ind_f1);
            f2_dir = f2_inv(ind_f1);
            
            figure;
            subplot(2,2,1)
            hist(f1_dir(find(f1_dir>0.03)))
            xlabel('F1')
            subplot(2,2,2)
            hist(f2_dir(find(f1_dir>0.03)))
            xlabel('F2')
            subplot(2,2,3)
            scatter(f1_dir(find(f1_dir>0.03)),f2_dir(find(f1_dir>0.03)))
            xlabel('F1')
            ylabel('F2')
            xlim([0 0.4])
            ylim([0 0.4])
            refline(1)
            subplot(2,2,4)
            f2overf1 = f2_dir./f1_dir;
            hist(f2overf1(find(f1_dir>0.03)))
            xlabel('F2/F1')
            suptitle([date ' ' mouse '- Phase reversal'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.pdf']),'-dpdf','-bestfit');
            save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']), 'f1', 'f2', 'f1_dir','f2_dir','f2overf1', 'resp_ind_phase')

            %% compare with adapt index
            
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ad_run_str], [date '_' mouse '_' ad_run_str '_dfofData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ad_run_str], [date '_' mouse '_' ad_run_str '_stimData.mat']))          
            resp_ind = find(sum(h1_ori,[2 3]));
            resp_pr = intersect(resp_ind, find(f1_dir>0.03));
            
            if ~exist('norm_resp_pref')
                norm_resp_pref = norm_dfof_stim_pref;
            end

            figure;
            scatter(f2overf1(resp_pr),norm_resp_pref(resp_pr),[],sfs(pref_sf(resp_pr)))
            xlim([0 1])
            ylim([-0.5 2])
            xlabel('F2/F1')
            ylabel('Norm resp')
            colorbar
            