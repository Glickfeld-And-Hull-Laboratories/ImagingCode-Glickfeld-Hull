%% get path names
close all;clear all;clc;

ds = 'F0F1F2_ExptList';
eval(ds)
nexp = length(expt);
rc = behavConstsAV;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

iexp = 15;
             %%
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        area = expt(iexp).img_loc{1};
        ImgFolder = expt(iexp).prFolder;
        if isfield(expt,'cofolder')
            coFolder = expt(iexp).coFolder;
        else
            coFolder = expt(iexp).dirFolder;
        end
        time = expt(iexp).prTime;
        nrun = length(ImgFolder);
        frameRateHz = params.frameRate;

        run_str = catRunName(cell2mat(ImgFolder), nrun);
        co_run_str = catRunName(cell2mat(coFolder), nrun);

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

                %% register to cross-ori experiment

                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_reg_shifts.mat']))
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
                title('Cross-ori run avg')
                sz = size(data_avg);
                rgb = zeros(sz(1),sz(2),3);
                rgb(:,:,1) = data_reg_avg./max(data_reg_avg(:));
                rgb(:,:,2) = data_avg./max(data_avg(:));
                subplot(2,2,3);
                image(rgb)
                title('Overlay')

                print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

                %% use cross-ori mask to get TCs

                fprintf(['Loading masks from cross-ori runs: ' cell2mat(coFolder) '\n'])

                % loads 'mask_cell', 'mask_np'
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_mask_cell.mat']))
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
                            [h_dir(:,iPhase,iDir,i), p_dir(:,iPhase,iDir,i)] = ttest2(squeeze(mean(data_dfof_cycavg(resp_win{i},:,ind),1)), squeeze(mean(data_dfof_cycavg(base_win{i},:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nDir-1));
                        end
                    end
                end
            end

            [max_val, max_ind] = max(max(max(phasedir_resp_avg(:,:,:,1),[],4),[],2),[],3);

            resp_ind_phase = find(sum(sum(sum(h_dir,2),3),4));
            figure;
            movegui('center')
            [n n2] = subplotn(length(resp_ind_phase));
            for iCell = 1:length(resp_ind_phase)
                subplot(n,n2,iCell)
                plot(tt, squeeze(data_dfof_phasedir(:,resp_ind_phase(iCell),:,max_ind(resp_ind_phase(iCell)))))
                hold on
                vline(tt(resp_win{1}))
                vline(tt(resp_win{2}))
                title([num2str(resp_ind_phase(iCell))])
            end
            suptitle([date ' ' mouse '- Phase reversal'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgRespByPhase_prefDir.pdf']),'-dpdf','-fillpage');

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
            f1 = zeros(1,nCells);
            f2 = zeros(1,nCells);
            for iCell = 1:nCells
                cyc = squeeze(data_dfof_phasedir(phaseCyc./2+1:end,iCell,:,max_ind(iCell)))';
                [f1(1,iCell),f2(1,iCell),f1ang,projectedf1,f1mat,f2mat] = compcontrastrevf1f2(cyc);
            end

            figure;
            subplot(2,2,1)
            hist(f1(resp_ind_phase))
            xlabel('F1')
            subplot(2,2,2)
            hist(f2(resp_ind_phase))
            xlabel('F2')
            subplot(2,2,3)
            scatter(f1(resp_ind_phase),f2(resp_ind_phase))
            xlabel('F1')
            ylabel('F2')
            xlim([0 0.4])
            ylim([0 0.4])
            refline(1)
            subplot(2,2,4)
            f2overf1 = f2./f1;
            hist(f2overf1(resp_ind_phase))
            xlabel('F2/F1')
            suptitle([date ' ' mouse '- Phase reversal'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.pdf']),'-dpdf','-bestfit');
            save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']), 'f1', 'f2', 'f2overf1', 'resp_ind_phase')

            %% compare with suppression index
            if isfield(expt,'cofolder')
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_respData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_dirAnalysis.mat']))          
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_dataStim.mat']))          
            
            resp_ind_only = intersect(resp_ind,resp_ind_phase);
            plaid_resp = mean(resp_cell{end,end,1},2);
            mask_resp =  mean(resp_cell{end,1,1},2);
            test_resp =  mean(resp_cell{1,end,1},2);
            plaid_resp(find(plaid_resp<0)) = 0;
            mask_resp(find(mask_resp<0)) = 0;
            test_resp(find(test_resp<0)) = 0;
            plaidSI = (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp);
            testPI = abs((test_resp-mask_resp) ./ (mask_resp+test_resp));
            figure; 
            movegui('center')
            subplot(2,2,1)
            scatter(plaidSI(resp_ind_only),f2overf1(resp_ind_only))
            xlim([-1 1])
            ylim([0 2.5])
            xlabel('Suppression index')
            ylabel('F2/F1')
            subplot(2,2,2)
            scatter(plaidSI(resp_ind_only),testPI(resp_ind_only))
            xlim([-1 1])
            ylim([0 1])
            xlabel('Suppression index')
            ylabel('Stim Pref index')
            subplot(2,2,3)
            scatter(testPI(resp_ind_only),f2overf1(resp_ind_only))
            xlim([0 1])
            ylim([0 2.5])
            xlabel('Stim Pref index')
            ylabel('F2/F1')
            subplot(2,2,4)
            scatter(Rc(resp_ind_only),f2overf1(resp_ind_only))
            xlim([-1 10])
            ylim([0 2.5])
            xlabel('Rc')
            ylabel('F2/F1')
            suptitle([date ' ' mouse '- Phase reversal'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SIvsF2-F1.pdf']),'-dpdf','-bestfit');

            complex_ind = intersect(resp_ind_only,find(f2overf1>1));
            simple_ind = intersect(resp_ind_only,find(f2overf1<0.5));

            figure;
            movegui('center')
            [n n2] = subplotn(length(complex_ind));
            for iC = 1:length(complex_ind)
                subplot(n,n2,iC)
                errorbar(stimDirs,avg_resp_dir(complex_ind(iC),:,1,1), avg_resp_dir(complex_ind(iC),:,1,2),'o')
                hold on
                errorbar(stimDirs,avg_resp_dir(complex_ind(iC),:,2,1), avg_resp_dir(complex_ind(iC),:,2,2),'o')
                title(['Rc ' num2str(Rc(complex_ind(iC))) '; Rp ' num2str(Rp(complex_ind(iC)))])
            end
            suptitle([date ' ' mouse '- Complex Cells'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_complexTuning.pdf']),'-dpdf','-bestfit');


            figure;
            movegui('center')
            [n n2] = subplotn(length(simple_ind));
            for iC = 1:length(simple_ind)
                subplot(n,n2,iC)
                errorbar(stimDirs,avg_resp_dir(simple_ind(iC),:,1,1), avg_resp_dir(simple_ind(iC),:,1,2),'o')
                hold on
                errorbar(stimDirs,avg_resp_dir(simple_ind(iC),:,2,1), avg_resp_dir(simple_ind(iC),:,2,2),'o')
                title(['Rc ' num2str(Rc(simple_ind(iC))) '; Rp ' num2str(Rp(simple_ind(iC)))])
            end
            suptitle([date ' ' mouse '- Simple Cells'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_simpleTuning.pdf']),'-dpdf','-bestfit');
            end
