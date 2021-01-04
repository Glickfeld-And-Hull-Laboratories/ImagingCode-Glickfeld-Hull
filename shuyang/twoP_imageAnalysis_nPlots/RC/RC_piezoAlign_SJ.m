clear all
close all

analysis_out = 'Z:\2P_analysis\';
bdata_source = 'Z:\Data\behavior\RC\';
CRP_expt_list_mike;
pairings = Load_CRP_List();


for  j = 1:size(pairings,2)
id = pairings{1,j};
exp_subset = pairings{2,j};
nexp = size(expt(id).date,1);
   for iexp = exp_subset
        mouse = expt(id).mouse(iexp,:);
        if find(mouse == ' ')
            ind = find(mouse == ' ');
            mouse(ind) = [];
        end
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        
        cd(fullfile('H:\home\mike\Data\2P Imaging', [date '_img' mouse],['img' mouse]))
        load(['img' mouse '_000_' run '.mat']);
        fn_piezo = fopen(['img' mouse '_000_' run '.ephys']);
        piezo_data = fread(fn_piezo,'single');
        first_frame = find(piezo_data==1,1,'first');
        nf = info.config.frames;
        nf = 85651;
        last_frame = find(piezo_data==nf,1,'first')+5;
        piezo_data_temp = piezo_data(first_frame:last_frame);
        piezo_data_temp = reshape(piezo_data_temp,[2 size(piezo_data_temp,1)./2]);
        piezo_frames = zeros(nf,1);
        for iframe = 1:nf
            ind = find(piezo_data_temp(1,:)==iframe);
            piezo_frames(iframe,:) = mean(piezo_data_temp(2,ind),2);
        end
        
        load(fullfile(analysis_out,img_fn, [img_fn '_input.mat']));
        load(fullfile(analysis_out,img_fn, [img_fn '_targetAlign.mat']));
        
        if id == 4
            rewDelay_frames =  round(1.1.*frameRateHz);
        else
            rewDelay_frames =  round(0.6.*frameRateHz);
        end
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        nTrials = length(cTargetOn);
        targetAlign_piezo = nan(prewin_frames+postwin_frames,nTrials);
        for itrial = 1:nTrials
            if cTargetOn(itrial)+postwin_frames-1 <= nf 
                targetAlign_piezo(:,itrial) = piezo_frames(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
            end
        end
        
        ind_nan = find(isnan(targetAlign_piezo(1,:)));
        
        
        targetAlign_piezo = abs(targetAlign_piezo);
        load(fullfile(analysis_out,img_fn, [img_fn '_cueAlignLick.mat']))
        figure;
        subplot(2,1,1)
        shadedErrorBar(tt,nanmean(targetAlign_piezo(:,ind_rew),2),nanstd(targetAlign_piezo(:,ind_rew),[],2)./sqrt(length(ind_rew)),'k');
        hold on
        if length(ind_omit>5)
            shadedErrorBar(tt,nanmean(targetAlign_piezo(:,ind_omit),2),nanstd(targetAlign_piezo(:,ind_omit),[],2)./sqrt(length(ind_omit)),'r');
        end
        if length(ind_unexp>5)
            shadedErrorBar(tt,nanmean(targetAlign_piezo(:,ind_unexp),2),nanstd(targetAlign_piezo(:,ind_unexp),[],2)./sqrt(length(ind_unexp)),'r');
        end
        ylabel('Piezo voltage')
        xlabel('Time from cue (ms)')
        subplot(2,1,2)
        shadedErrorBar(tt,nanmean(lickCueAlign(:,ind_rew),2),nanstd(lickCueAlign(:,ind_rew),[],2)./sqrt(length(ind_rew)),'k');
        hold on
        if length(ind_omit>5)
            shadedErrorBar(tt,nanmean(lickCueAlign(:,ind_omit),2),nanstd(lickCueAlign(:,ind_omit),[],2)./sqrt(length(ind_omit)),'r');
        end
        if length(ind_unexp>5)
            shadedErrorBar(tt,nanmean(lickCueAlign(:,ind_unexp),2),nanstd(lickCueAlign(:,ind_unexp),[],2)./sqrt(length(ind_unexp)),'r');
        end
        ylabel('Lick rate')
        xlabel('Time from cue (ms)')
        suptitle([date ' ' mouse '- Reward (black), Omit (red)'])
        savefig(fullfile(analysis_out,img_fn, [img_fn '_avgTrialPiezo_abs.fig']))
        
        load(fullfile(analysis_out,img_fn, [img_fn '_cueAlignLick.mat']))
        preRew_lickSearchRange = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
        postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:prewin_frames+rewDelay_frames+rewDelay_frames;
        preRew_piezoAmp = mean(targetAlign_piezo(preRew_lickSearchRange,:),1);
        postRew_piezoAmp = mean(targetAlign_piezo(postRew_lickSearchRange,:),1);
        
        figure;
        subplot(1,2,1)
        scatter(preRew_lickBurstHz(:,ind_rew), preRew_piezoAmp(:,ind_rew),'ok')
        hold on
        scatter(preRew_lickBurstHz(:,ind_omit), preRew_piezoAmp(:,ind_omit),'or')
        xlabel('Lick rate')
        ylabel('Piezo voltage')
        title('Pre reward')
        xlim([0 10])
        ylim([-0.5 0.5])
        axis square
        subplot(1,2,2)
        scatter(postRew_lickBurstHz(:,ind_rew), postRew_piezoAmp(:,ind_rew),'ok')
        hold on
        scatter(postRew_lickBurstHz(:,ind_omit), postRew_piezoAmp(:,ind_omit),'or')
        xlabel('Lick rate')
        ylabel('Piezo voltage')
        title('Pre reward')
        xlim([0 10])
        ylim([-0.5 0.5])
        axis square
        suptitle([date ' ' mouse '- Reward (black), Omit (red)'])
        savefig(fullfile(analysis_out,img_fn, [img_fn '_LickvsPiezo_abs.fig']))
        
        [sortPiezoAmp sortPiezoAmp_ind] = sort(preRew_piezoAmp(:,ind_rew),'ascend');
        nnan = sum(isnan(preRew_piezoAmp(:,ind_rew)));
        ind_low25piezo_prerew = sortPiezoAmp_ind(1:floor(length(ind_rew)/4));
        ind_high25piezo_prerew = sortPiezoAmp_ind(length(ind_rew)-floor(length(ind_rew)/4)+1:end-nnan);
        HL_piezo.low25_prerew = nanmean(preRew_piezoAmp(:,ind_rew(ind_low25piezo_prerew)),2);
        HL_piezo.high25_prerew = nanmean(preRew_piezoAmp(:,ind_rew(ind_high25piezo_prerew)),2);
        ind_low10piezo_prerew = sortPiezoAmp_ind(1:floor(length(ind_rew)/10));
        ind_high10piezo_prerew = sortPiezoAmp_ind(length(ind_rew)-floor(length(ind_rew)/10)+1:end-nnan);
        HL_piezo.low10_prerew = nanmean(preRew_piezoAmp(:,ind_rew(ind_low10piezo_prerew)),2);
        HL_piezo.high10_prerew = nanmean(preRew_piezoAmp(:,ind_rew(ind_high10piezo_prerew)),2);

        [sortPiezoAmp sortPiezoAmp_ind] = sort(postRew_piezoAmp(:,ind_rew),'ascend');
        nnan = sum(isnan(postRew_piezoAmp(:,ind_rew)));
        ind_low25piezo_postrew = sortPiezoAmp_ind(1:floor(length(ind_rew)/4));
        ind_high25piezo_postrew = sortPiezoAmp_ind(length(ind_rew)-floor(length(ind_rew)/4)+1:end-nnan);
        HL_piezo.low25_postrew = nanmean(postRew_piezoAmp(:,ind_rew(ind_low25piezo_postrew)),2);
        HL_piezo.high25_postrew = nanmean(postRew_piezoAmp(:,ind_rew(ind_high25piezo_postrew)),2);
        ind_low10piezo_postrew = sortPiezoAmp_ind(1:floor(length(ind_rew)/10));
        ind_high10piezo_postrew = sortPiezoAmp_ind(length(ind_rew)-floor(length(ind_rew)/10)+1:end-nnan);
        HL_piezo.low10_postrew = nanmean(postRew_piezoAmp(:,ind_rew(ind_low10piezo_postrew)),2);
        HL_piezo.high10_postrew = nanmean(postRew_piezoAmp(:,ind_rew(ind_high10piezo_postrew)),2);
        if length(ind_omit)>5
            [sortPiezoAmp sortPiezoAmp_ind] = sort(preRew_piezoAmp(:,ind_omit),'ascend');
            nnan = sum(isnan(preRew_piezoAmp(:,ind_omit)));
            ind_low25piezo_preomit = sortPiezoAmp_ind(1:floor(length(ind_omit)/4));
            ind_high25piezo_preomit = sortPiezoAmp_ind(length(ind_omit)-floor(length(ind_omit)/4)+1:end-nnan);
            HL_piezo.low25_preomit = nanmean(preRew_piezoAmp(:,ind_omit(ind_low25piezo_preomit)),2);
            HL_piezo.high25_preomit = nanmean(preRew_piezoAmp(:,ind_omit(ind_high25piezo_preomit)),2);
            ind_low10piezo_preomit = sortPiezoAmp_ind(1:floor(length(ind_omit)/10));
            ind_high10piezo_preomit = sortPiezoAmp_ind(length(ind_omit)-floor(length(ind_omit)/10)+1:end-nnan);
            HL_piezo.low10_preomit = nanmean(preRew_piezoAmp(:,ind_omit(ind_low10piezo_preomit)),2);
            HL_piezo.high10_preomit = nanmean(preRew_piezoAmp(:,ind_omit(ind_high10piezo_preomit)),2);

            [sortPiezoAmp sortPiezoAmp_ind] = sort(postRew_piezoAmp(:,ind_omit),'ascend');
            nnan = sum(isnan(postRew_piezoAmp(:,ind_omit)));
            ind_low25piezo_postomit = sortPiezoAmp_ind(1:floor(length(ind_omit)/4));
            ind_high25piezo_postomit = sortPiezoAmp_ind(length(ind_omit)-floor(length(ind_omit)/4)+1:end-nnan);
            HL_piezo.low25_postomit = nanmean(postRew_piezoAmp(:,ind_omit(ind_low25piezo_postomit)),2);
            HL_piezo.high25_postomit = nanmean(postRew_piezoAmp(:,ind_omit(ind_high25piezo_postomit)),2);
            ind_low10piezo_postomit = sortPiezoAmp_ind(1:floor(length(ind_omit)/10));
            ind_high10piezo_postomit = sortPiezoAmp_ind(length(ind_omit)-floor(length(ind_omit)/10)+1:end-nnan);
            HL_piezo.low10_postomit = nanmean(postRew_piezoAmp(:,ind_omit(ind_low10piezo_postomit)),2);
            HL_piezo.high10_postomit = nanmean(postRew_piezoAmp(:,ind_omit(ind_high10piezo_postomit)),2);
        else
            ind_low25piezo_preomit = [];
            ind_high25piezo_preomit = [];
            ind_low25piezo_postomit = [];
            ind_high25piezo_postomit =[];
            HL_piezo.low25_preomit = [];
            HL_piezo.high25_preomit = [];
            HL_piezo.low25_postomit = [];
            HL_piezo.high25_postomit = [];
            ind_low10piezo_preomit = [];
            ind_high10piezo_preomit = [];
            ind_low10piezo_postomit = [];
            ind_high10piezo_postomit =[];
            HL_piezo.low10_preomit = [];
            HL_piezo.high10_preomit = [];
            HL_piezo.low10_postomit = [];
            HL_piezo.high10_postomit = [];
        end
        nIC = size(targetAlign_events,2);
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_low25piezo_prerew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_low25piezo_prerew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
        hold on
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_high25piezo_prerew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_high25piezo_prerew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 inf])
        title(['Pre- Rew: ' num2str(chop(HL_piezo.low25_prerew,2)) ' vs ' num2str(chop(HL_piezo.high25_prerew,2)) ' V'])
        subplot(2,2,3)
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_low25piezo_postrew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_low25piezo_postrew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
        hold on
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_high25piezo_postrew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_high25piezo_postrew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 inf])
        title(['Post- Rew: ' num2str(chop(HL_piezo.low25_postrew,2)) ' vs ' num2str(chop(HL_piezo.high25_postrew,2)) ' V'])
        if length(ind_omit)>5
            subplot(2,2,2)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_low25piezo_preomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_low25piezo_preomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
            hold on
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_high25piezo_preomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_high25piezo_preomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 inf])
            title(['Pre- Omit: ' num2str(chop(HL_piezo.low25_preomit,2)) ' vs ' num2str(chop(HL_piezo.high25_preomit,2)) ' V'])
            subplot(2,2,4)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_low25piezo_postomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_low25piezo_postomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
            hold on
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_high25piezo_postomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_high25piezo_postomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 inf])
            title(['Post- Omit: ' num2str(chop(HL_piezo.low25_postomit,2)) ' vs ' num2str(chop(HL_piezo.high25_postomit,2)) ' V'])
        end
        suptitle([mouse ' ' date '- lick bursts by rate: low25 (blue) & high25 (black)'])
        savefig(fullfile(analysis_out,img_fn, [img_fn '_cueAlignSpiking_byPiezoAmp25_abs.fig'])) 
        
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_low10piezo_prerew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_low10piezo_prerew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
        hold on
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_high10piezo_prerew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_high10piezo_prerew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 inf])
        title(['Pre- Rew: ' num2str(chop(HL_piezo.low10_prerew,2)) ' vs ' num2str(chop(HL_piezo.high10_prerew,2)) ' V'])
        subplot(2,2,3)
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_low10piezo_postrew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_low10piezo_postrew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
        hold on
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_high10piezo_postrew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_high10piezo_postrew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 inf])
        title(['Post- Rew: ' num2str(chop(HL_piezo.low10_postrew,2)) ' vs ' num2str(chop(HL_piezo.high10_postrew,2)) ' V'])
        if length(ind_omit)>5
            subplot(2,2,2)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_low10piezo_preomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_low10piezo_preomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
            hold on
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_high10piezo_preomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_high10piezo_preomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 inf])
            title(['Pre- Omit: ' num2str(chop(HL_piezo.low10_preomit,2)) ' vs ' num2str(chop(HL_piezo.high10_preomit,2)) ' V'])
            subplot(2,2,4)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_low10piezo_postomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_low10piezo_postomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
            hold on
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_high10piezo_postomit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_high10piezo_postomit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 inf])
            title(['Post- Omit: ' num2str(chop(HL_piezo.low10_postomit,2)) ' vs ' num2str(chop(HL_piezo.high10_postomit,2)) ' V'])
        end
        suptitle([mouse ' ' date '- lick bursts by rate: low10 (blue) & high10 (black)'])
        savefig(fullfile(analysis_out,img_fn, [img_fn '_cueAlignSpiking_byPiezoAmp10_abs.fig'])) 
        
        piezo_base_std = nanstd(reshape(targetAlign_piezo(1:prewin_frames,:),[prewin_frames.*nTrials 1]),[],1);
        piezo_base_avg = nanmean(reshape(targetAlign_piezo(1:prewin_frames,:),[prewin_frames.*nTrials 1]),1);
        targetAlign_piezo_thresh1 = zeros(size(targetAlign_piezo));
        targetAlign_piezo_thresh1ALL = zeros(size(targetAlign_piezo));
        targetAlign_piezo_thresh2ALL = zeros(size(targetAlign_piezo));
        targetAlign_piezo_thresh2 = zeros(size(targetAlign_piezo));
        targetAlign_piezo_thresh3 = zeros(size(targetAlign_piezo));
        piezoSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
        piezoStart = nan(1,nTrials);
        for itrial = 1:nTrials
            targetAlign_piezo_thresh1(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*1) & targetAlign_piezo(:,itrial)<(piezo_base_avg + piezo_base_std.*2)),itrial) = 1;
            targetAlign_piezo_thresh2(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*2) & targetAlign_piezo(:,itrial)<(piezo_base_avg + piezo_base_std.*3)),itrial) = 1;
            targetAlign_piezo_thresh3(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*3)),itrial) = 1;
            targetAlign_piezo_thresh1ALL(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*1)),itrial) = 1;
            targetAlign_piezo_thresh2ALL(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*2)),itrial) = 1;
            if find(targetAlign_piezo_thresh3(piezoSearchRange,itrial),1,'first')
                piezoStart(:,itrial) = find(targetAlign_piezo_thresh3(piezoSearchRange,itrial),1,'first');
            end
        end
        
        [sortPiezoStart sortPiezoStart_ind] = sort(piezoStart(:,ind_rew),'ascend');
        nnan = sum(isnan(piezoStart(:,ind_rew)));
        ind_earlypiezo_rew = sortPiezoStart_ind(1:floor((length(sortPiezoStart_ind)-nnan)/4));
        ind_latepiezo_rew = sortPiezoStart_ind((length(sortPiezoStart_ind)-nnan)-floor((length(sortPiezoStart_ind)-nnan)/4)+1:end-nnan);
        ind_allearlypiezo25_rew = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/4));
        ind_alllatepiezo25_rew = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/4)+1:end);
        ind_allearlypiezo10_rew = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/10));
        ind_alllatepiezo10_rew = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/10)+1:end);
        HL_piezo.early_rew = nanmean(piezoStart(:,ind_rew(ind_earlypiezo_rew)),2);
        HL_piezo.late_rew = nanmean(piezoStart(:,ind_rew(ind_latepiezo_rew)),2);
        if length(ind_omit)>5
            [sortPiezoStart sortPiezoStart_ind] = sort(piezoStart(:,ind_omit),'ascend');
            nnan = sum(isnan(piezoStart(:,ind_omit)));
            ind_earlypiezo_omit = sortPiezoStart_ind(1:floor((length(sortPiezoStart_ind)-nnan)/4));
            ind_latepiezo_omit = sortPiezoStart_ind((length(sortPiezoStart_ind)-nnan)-floor((length(sortPiezoStart_ind)-nnan)/4)+1:end-nnan);
            ind_allearlypiezo25_omit = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/4));
            ind_alllatepiezo25_omit = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/4)+1:end);
            ind_allearlypiezo10_omit = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/10));
            ind_alllatepiezo10_omit = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/10)+1:end);
            HL_piezo.early_omit = nanmean(piezoStart(:,ind_omit(ind_earlypiezo_omit)),2);
            HL_piezo.late_omit = nanmean(piezoStart(:,ind_omit(ind_latepiezo_omit)),2);
        else
            ind_earlypiezo_omit = [];
            ind_latepiezo_omit = [];
            ind_allearlypiezo10_omit = [];
            ind_alllatepiezo10_omit = [];
            ind_allearlypiezo25_omit = [];
            ind_alllatepiezo25_omit = [];
            HL_piezo.late_omit = [];
            HL_piezo.early_omit = [];
        end
        
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_earlypiezo_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_earlypiezo_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
        hold on
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_latepiezo_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_latepiezo_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        scatter(piezoStart(:,ind_rew(ind_earlypiezo_rew)).*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew)),'ob')
        scatter(piezoStart(:,ind_rew(ind_latepiezo_rew)).*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew)),'ok')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 inf])
        title(['Rew: ' num2str(chop(HL_piezo.early_rew.*(1000./frameRateHz),2)) ' vs ' num2str(chop(HL_piezo.late_rew.*(1000./frameRateHz),2)) ' ms'])
        subplot(2,2,3)
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_allearlypiezo10_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_allearlypiezo10_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
        hold on
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_alllatepiezo10_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_alllatepiezo10_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 inf])
        title(['All- Rew'])
        if length(ind_omit)>5
            subplot(2,2,2)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_earlypiezo_omit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_earlypiezo_omit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
            hold on
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_latepiezo_omit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_latepiezo_omit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            scatter(piezoStart(:,ind_omit(ind_earlypiezo_omit)).*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_omit)),'ob')
            scatter(piezoStart(:,ind_omit(ind_latepiezo_omit)).*(1000./frameRateHz),zeros(1,length(ind_latepiezo_omit)),'ok')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 inf])
            title(['All- Omit: ' num2str(chop(HL_piezo.early_omit.*(1000./frameRateHz),2)) ' vs ' num2str(chop(HL_piezo.late_omit.*(1000./frameRateHz),2)) ' ms'])
            subplot(2,2,4)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_allearlypiezo10_omit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_allearlypiezo10_omit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
            hold on
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_alllatepiezo10_omit)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_alllatepiezo10_omit)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 inf])
            title(['All- Omit'])
        end
        suptitle([mouse ' ' date '- lick bursts by latency: early (blue) & late (black)'])
        savefig(fullfile(analysis_out,img_fn, [img_fn '_cueAlignSpiking_byPiezoLatency_abs.fig'])) 
        
        save(fullfile(analysis_out,img_fn, [img_fn '_cueAlignPiezo.mat']), 'targetAlign_piezo','targetAlign_piezo_thresh1','targetAlign_piezo_thresh1ALL','targetAlign_piezo_thresh2','targetAlign_piezo_thresh2ALL','targetAlign_piezo_thresh3', 'preRew_piezoAmp', 'postRew_piezoAmp', 'ind_low25piezo_prerew', 'ind_high25piezo_prerew','ind_low25piezo_postrew', 'ind_high25piezo_postrew', 'ind_low25piezo_preomit', 'ind_high25piezo_preomit','ind_low25piezo_postomit', 'ind_high25piezo_postomit','ind_low10piezo_prerew', 'ind_high10piezo_prerew','ind_low10piezo_postrew', 'ind_high10piezo_postrew', 'ind_low10piezo_preomit', 'ind_high10piezo_preomit','ind_low10piezo_postomit', 'ind_high10piezo_postomit', 'ind_earlypiezo_rew', 'ind_latepiezo_rew','ind_earlypiezo_omit', 'ind_latepiezo_omit', 'ind_allearlypiezo25_rew', 'ind_alllatepiezo25_rew','ind_allearlypiezo25_omit', 'ind_alllatepiezo25_omit', 'ind_allearlypiezo10_rew', 'ind_alllatepiezo10_rew','ind_allearlypiezo10_omit', 'ind_alllatepiezo10_omit', 'HL_piezo')
    end
end
              
        