homeDir = 'Z:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis';
xlsfile = fullfile(homeDir, 'LED-full.xlsx');
[Gdata, Gtext, Graw] = xlsread(xlsfile,23);
IDs =Gdata(:,1);
Dir = Gtext(2:end,5);
region = Gtext(2:end,2);
power = Gdata(:,3);
background = Gtext(2:end,6);
Background = Gtext(2:end,8);
virus = Gtext(2:end,7);
spread = Gdata(:,9);
Thresh = NaN(length(IDs),2); %  first collumn is the control, second is led
FA =  NaN(length(IDs),2);
FA_CI = NaN(length(IDs),2,2); % confidence interval for each mice
Thresh_CI = NaN(length(IDs),2,2); % confidence interval for each mice
slope = NaN(length(IDs),2);
Lapse = NaN(length(IDs),2); % store lapse rate
Fidget = NaN(length(IDs),2);
FA_trial = NaN(length(IDs),2);
dprime = NaN(length(IDs),12,2);
criterion = NaN(length(IDs),12,2);
Oriens = NaN(length(IDs),12,2);
RT = NaN(length(IDs),12,2); % store the reaction time

LEDpre_thresh = NaN(length(IDs),4); % all combined: LED/Ctrl; LED/LED; Ctrl/Ctrl; Ctrl/LED
LEDpre_FA = NaN(length(IDs),4); % previouse led trial, followed by control/or followed by led
LEDpre_hit = NaN(length(IDs),4);% for interpolated 22.5 deg
LEDpre_d =  NaN(length(IDs),4);
LEDpre_c = NaN(length(IDs),4);
LEDpre_slope = NaN(length(IDs),4);
LEDpre_lapse = NaN(length(IDs),4);


FA_bintrial = NaN(length(IDs),4,2); % only has four bins
Thresh_bintrial = NaN(length(IDs),4,2);
bootStats = {}; % for each individual mice

Hit_22 = NaN(length(IDs),2);
d_22 = NaN(length(IDs),2);
c_22 = NaN(length(IDs),2);

for i_exp = 1:length(IDs)
    load(Dir{i_exp})
    
    % get the threshold for previouse LED, current trial is control (LED/Ctrl)
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.LED_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.LED_pre.all.HT_num + Output{1}.Outcome{1}.target.LED_pre.all.Miss_num;
    LEDpre_fit{i_exp,1} = weibullFitLG(Output{1}.Outcome{1}.target.LED_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,1) = LEDpre_fit{i_exp,1}.thresh;
    LEDpre_FA(i_exp,1) = Output{1}.Outcome{1}.FA.LED_pre.all.FA;
    LEDpre_hit(i_exp,1) = LEDpre_fit{i_exp,1}.modelFun(LEDpre_fit{i_exp,1}.coefEsts, 22.5);
    [LEDpre_d(i_exp,1), LEDpre_c(i_exp,1)] = dprime_simple(LEDpre_hit(i_exp,1),LEDpre_FA(i_exp,1));
    LEDpre_lapse(i_exp,1) = 1-Hits(end); 
    LEDpre_slope(i_exp,1) = LEDpre_fit{i_exp,1}.slope; 
    
    % get the threshold for previouse LED, current trial is LED(LED/LED)
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.Same_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.Same_pre.all.HT_num + Output{1}.Outcome{2}.target.Same_pre.all.Miss_num;
    LEDpre_fit{i_exp,2} = weibullFitLG(Output{1}.Outcome{2}.target.Same_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,2) = LEDpre_fit{i_exp,2}.thresh;
    LEDpre_FA(i_exp,2) = Output{1}.Outcome{2}.FA.Same_pre.all.FA;
    LEDpre_hit(i_exp,2) = LEDpre_fit{i_exp,2}.modelFun(LEDpre_fit{i_exp,2}.coefEsts, 22.5); 
    [LEDpre_d(i_exp,2), LEDpre_c(i_exp,2)] = dprime_simple(LEDpre_hit(i_exp,2),LEDpre_FA(i_exp,2));
    LEDpre_lapse(i_exp,2) = 1-Hits(end); 
    LEDpre_slope(i_exp,2) = LEDpre_fit{i_exp,2}.slope; 
    
    % get the threshold for previouse control, current trial is control (Ctrl/Ctrl)
    Hits = [];
    Hits = Output{1}.Outcome{1}.target.Same_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{1}.target.Same_pre.all.HT_num + Output{1}.Outcome{1}.target.Same_pre.all.Miss_num;
    LEDpre_fit{i_exp,3} = weibullFitLG(Output{1}.Outcome{1}.target.Same_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,3) = LEDpre_fit{i_exp,3}.thresh;
    LEDpre_FA(i_exp,3) = Output{1}.Outcome{1}.FA.Same_pre.all.FA;
    LEDpre_hit(i_exp,3) = LEDpre_fit{i_exp,3}.modelFun(LEDpre_fit{i_exp,3}.coefEsts, 22.5);
    [LEDpre_d(i_exp,3), LEDpre_c(i_exp,3)] = dprime_simple(LEDpre_hit(i_exp,3),LEDpre_FA(i_exp,3));
    LEDpre_lapse(i_exp,3) = 1-Hits(end); 
    LEDpre_slope(i_exp,3) = LEDpre_fit{i_exp,3}.slope; 
    
    % get the threshold for previouse Control, current trial is LED(Ctrl/LED)
    Hits = [];
    Hits = Output{1}.Outcome{2}.target.LED_pre.all.c_hit;
    trialAll = [];
    trialAll = Output{1}.Outcome{2}.target.LED_pre.all.HT_num + Output{1}.Outcome{2}.target.LED_pre.all.Miss_num;
    LEDpre_fit{i_exp,4} = weibullFitLG(Output{1}.Outcome{2}.target.LED_pre.Orien, Hits',1, 1, {'nTrials', trialAll'});
    LEDpre_thresh(i_exp,4) = LEDpre_fit{i_exp,4}.thresh;
    LEDpre_FA(i_exp,4) = Output{1}.Outcome{2}.FA.LED_pre.all.FA;
    LEDpre_hit(i_exp,4) = LEDpre_fit{i_exp,4}.modelFun(LEDpre_fit{i_exp,4}.coefEsts, 22.5);
    [LEDpre_d(i_exp,4), LEDpre_c(i_exp,4)] = dprime_simple(LEDpre_hit(i_exp,4),LEDpre_FA(i_exp,4));
    LEDpre_lapse(i_exp,4) = 1-Hits(end); 
    LEDpre_slope(i_exp,4) = LEDpre_fit{i_exp,4}.slope; 
        
    
    for i = 1:2
        Hits = [];
        trialAll = [];
        Hits = Output{i}.target.all.c_hit;
        trialAll = Output{i}.target.all.HT_num + Output{i}.target.all.Miss_num;
        fitSall{i_exp,i} = weibullFitLG(Output{i}.Infor.Orien, Hits',1, 1, {'nTrials', trialAll'});
        % add bootstrap for each mouse
        bootStats{i_exp,i}=BootstrapWeibullFit(trialAll',Hits',1000,Output{i}.Infor.Orien,1, 1);
        Thresh(i_exp,i) = fitSall{i_exp,i}.thresh;
        Thresh_CI(i_exp,i,:) = bootStats{i_exp,i}.ci95; 
        slope(i_exp,i) = fitSall{i_exp,i}.slope; 
        FA(i_exp,i) = Output{i}.FA.all.FA;
        FA_CI(i_exp,i,:) = Output{i}.FA.all.FA_confi; 
        Fidget(i_exp,i) = Output{i}.FA_percent.fidget;
        FA_trial(i_exp,i) = Output{i}.FA_percent.FA;
        dprime(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.sdt.all.dprime;
        criterion(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.sdt.all.criterion;
        Oriens(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.Infor.Orien;
        RT(i_exp,1:length(Output{i}.Infor.Orien),i) = Output{i}.target.all.RTonHit_mean(1,:);
        Lapse(i_exp,i) = 1-Hits(end);
        
        FA_bintrial(i_exp,1:4,i)  = Output{i}.bintrial.FA;
        Thresh_bintrial(i_exp,1:4,i)  =  Output{i}.bintrial.thresh;
        % interpolate the 22.5 deg for d' and c
        Hit_22(i_exp,i) = fitSall{i_exp,i}.modelFun(fitSall{i_exp,i}.coefEsts, 22.5); 
        [d_22(i_exp,i), c_22(i_exp,i)] = dprime_simple(Hit_22(i_exp,i),FA(i_exp,i));
%         d_22(i_exp,i) = mean(Output{i}.sdt.all.dprime(Output{i}.Infor.Orien>20 & Output{i}.Infor.Orien<30));
%         c_22(i_exp,i) = mean(Output{i}.sdt.all.criterion(Output{i}.Infor.Orien>20 & Output{i}.Infor.Orien<30));
    end
end
a=squeeze(Oriens(:,:,1));
cbin = Output{1}.bintrial.cbin.Hit; 
%% determine if significant different from control and LED for individual mouse
FA_sig = [];
Thresh_sig = []; 
for i_exp =1:length(IDs)
    
    FA_sig(i_exp,1) = max(FA_CI(i_exp,1,:))<min(FA_CI(i_exp,2,:))|| max(FA_CI(i_exp,2,:))<min(FA_CI(i_exp,1,:));
    Thresh_sig(i_exp,1) =  max(Thresh_CI(i_exp,1,:))<min(Thresh_CI(i_exp,2,:))|| max(Thresh_CI(i_exp,2,:))<min(Thresh_CI(i_exp,1,:));
end

%% for add on analysis
FA_RT = NaN(length(IDs),2);
for i_exp = 1:length(IDs)
    load(Dir{i_exp})
    for i = 1:2
        if iscell(Output{i}.FA.FA_RT)
            FA_RT(i_exp,i)= mean([Output{i}.FA.FA_RT{:}]);
        else
            FA_RT(i_exp,i)= mean(Output{i}.FA.FA_RT);
        end
    end
end

RT_easy = [];
RT_thresh = []; 
for i = 1:size(Oriens,1)
% for easy trial, find the nearest value to the easiest control orien
temp1 = [];
temp1 = squeeze(Oriens(i,:,1)); 
temp2 = [];
temp2 = squeeze(Oriens(i,:,2)); 
[easy_ori(i), idx1] = max(temp1); 
[c, idx2] = min(abs(temp2-easy_ori(i))); 
RT_easy (i,1) = RT(i,idx1,1); 
RT_easy (i,2) = RT(i,idx2,2); 

[c, idx3] = min(abs(temp1-Thresh(i,1))); 
thresh_ori(i) = temp1(idx3); 
[c, idx4] = min(abs(temp2-thresh_ori(i))); 
RT_thresh(i,1) = RT(i,idx3,1);
RT_thresh(i,2) = RT(i,idx4,2);

end 
%% only select mouse that has at least 20 trials for each bin/ori 

Mice_select = []; % for trial bin analysis
for i_exp = 1:length(IDs)
    load(Dir{i_exp})
    trial_c=min( Output{1}.bintrial.Trial_num); 
    trial_l=min( Output{2}.bintrial.Trial_num);
    if length(trial_c)>5
       Mice_select(i_exp) = trial_c(2)>=20; 
    else
       Mice_select(i_exp) = trial_c(1)>=20; 
    end 
   
end 
%%  summary of LED effects, delta change
Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
all_region = {'V1','LM','AL','PM'};
figure
subplot(2,3,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);
[h_t(i_region),p_t(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')

hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
 text(i_region,19,['n=' num2str(exp)])
ylim([-5 20])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
axis square

subplot(2,3,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);
[h_f(i_region),p_f(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])
ylim([-0.1 0.1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
axis square

subplot(2,3,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Lapse(idx,:);
[h_f(i_region),p_f(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])


end 
ylim([-0.05 0.15])
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Lapse rate(LED - Ctrl)')
axis square

subplot(2,3,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = slope(idx,:);
[h_f(i_region),p_f(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])


end 
ylim([-0.7 1.4])
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Slope(LED - Ctrl)')
axis square

subplot(2,3,5)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = d_22(idx,:);
[h_d(i_region),p_d(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.7,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta d prime-22.5 (LED - Ctrl)')
axis square

subplot(2,3,6)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = c_22(idx,:);
[h_c(i_region),p_c(i_region)]= ttest(temp(:,1),temp(:,2));
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)),std(diff(temp,1,2))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%text(i_region,0.7,['n=' num2str(exp)])
ylim([-1 1])

end 

xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta c-22.5 (LED - Ctrl)')
axis square

%% plot threshold over threshold
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);

scatter(temp(:,1),temp(:,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
hold on


end 
ylim([0 45])
xlim([0 45])
hold on
plot([0 45],[0 45],'k:')
set(gca,'TickDir','out')
ylabel('Thresh-LED')
xlabel('Thresh-Ctrl')
axis square
subplot(2,2,2)

for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);
scatter(temp(:,1),diff(temp,[],2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
hold on


end 
ylim([-5 20])
xlim([0 45])
set(gca,'TickDir','out')
ylabel('deltaThresh(LED-Ctrl)')
xlabel('Thresh-Ctrl')
axis square

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);

scatter(temp(:,1),temp(:,2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
hold on


end 
ylim([0 0.2])
xlim([0 0.2])
hold on
plot([0 0.2],[0 0.2],'k:')
set(gca,'TickDir','out')
ylabel('FA-LED')
xlabel('FA-Ctrl')
axis square
subplot(2,2,4)

for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);
scatter(temp(:,1),diff(temp,[],2),'MarkerEdgeColor',Colors(i_region,:),'SizeData',40)
hold on


end 
ylim([-0.1 0.1])
xlim([0 0.2])
set(gca,'TickDir','out')
ylabel('deltaFA(LED-Ctrl)')
xlabel('FA-Ctrl')
axis square
%% plot within each area, LED effects seperate by suppression method
back_color = [1 0.6 0.6;0 0 0; 0.5 0.5 0.5; 0.6 0.8 0.8];  % black PV-ChR2, gray PV-Chronos, red:GAD-ChR2; green:vgat-ChR2
B_unique = unique(Background); 

figure
for i_region = 1:4
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region));
    for i_b =1:4
        idx2 = strcmp(Background,B_unique(i_b));
        idx3 = idx&idx2;
        temp = Thresh(idx3,:);
        if ~isempty(temp)
            scatter(temp(:,1),temp(:,2),'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)
            hold on
        end
    end
    text(1,40,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 45])
    xlim([0 45])
    hold on
    plot([0 45],[0 45],'k:')
    set(gca,'TickDir','out','XTick',0:15:45, 'YTick',0:15:45)
    ylabel('Thresh-LED')
    xlabel('Thresh-Ctrl')
    axis square
    
    
end

figure

for i_region = 1:4
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region));
    for i_b =1:4
        idx2 = strcmp(Background,B_unique(i_b));
        idx3 = idx&idx2;
        temp = FA(idx3,:);
        if ~isempty(temp)
            scatter(temp(:,1),temp(:,2),'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)
            hold on
        end
    end
    text(0.02,0.18,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 0.2])
    xlim([0 0.2])
    hold on
    plot([0 0.2],[0 0.2],'k:')
    set(gca,'TickDir','out','XTick',0:0.1:0.2, 'YTick',0:0.1:0.2)
    ylabel('FA-LED')
    xlabel('FA-Ctrl')
    axis square
    
    
end

%% plot delta change over light power

figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region))& Thresh_sig;
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)

hold on

idx = strcmp(region,all_region(i_region))& ~Thresh_sig;
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',40)



end 
ylim([-5 20])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
xlabel('light power (mW)')
axis square

subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region))& FA_sig;
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)

hold on
idx = strcmp(region,all_region(i_region))& ~FA_sig;
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',40)


end 
ylim([-0.1 0.1])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
xlabel('light power (mW)')
axis square

subplot(2,2,3)
for i_b = 1:4
idx = strcmp(Background,B_unique(i_b))& Thresh_sig;
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on

idx = strcmp(Background,B_unique(i_b))& ~Thresh_sig;
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',[1 1 1],'SizeData',40)


end 
ylim([-5 20])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
xlabel('light power (mW)')
axis square

subplot(2,2,4)
for i_b = 1:4
idx = strcmp(Background,B_unique(i_b))& FA_sig;
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on
idx = strcmp(Background,B_unique(i_b))& ~FA_sig;
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',[1 1 1],'SizeData',40)


end 
ylim([-0.1 0.1])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
xlabel('light power (mW)')
axis square



% subplot(2,2,3)
% for i_region = 1:4
% idx = strcmp(region,all_region(i_region));
% temp = Lapse(idx,:);
% Temp = diff(temp,[],2);
% 
% scatter(power(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% 
% hold on
% 
% 
% end 
% ylim([-0.06 0.12])
% xlim([0 1])
% hline(0,'k:')
% set(gca,'TickDir','out','YTick',[-0.06:0.06:0.12])
% ylabel('\Delta Lapse (LED - Ctrl)')
% xlabel('light power (mW)')
% axis square
%% get the significant change in lapse rate
 Lapse_new = [];
 Lapse_CI = []; 
 Lapse_sig = []; 
for i_exp = 1:length(IDs)
    load(Dir{i_exp})
    for i = 1:2
       [ Lapse_new(i_exp,i), Lapse_CI(i_exp,i,:)] = binofit(Output{i}.target.all.Miss_num(end),Output{i}.target.all.Miss_num(end)+Output{i}.target.all.HT_num (end));  
    
    end
     Lapse_sig(i_exp,1) = max(Lapse_CI(i_exp,1,:))<min(Lapse_CI(i_exp,2,:))|| max(Lapse_CI(i_exp,2,:))<min(Lapse_CI(i_exp,1,:));
    
end


%% plot the delta change in thresh and FA over lapse rate changes

figure
D_lapse = diff(Lapse,[],2); 
subplot(1,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region))& Lapse_sig;
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(D_lapse(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)

hold on
idx = strcmp(region,all_region(i_region))& ~Lapse_sig;
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(D_lapse(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',40)

end 
% get the correlation coeficient 
[R,P] = corrcoef(diff(Thresh,[],2), diff(Lapse,[],2) ); 
ylim([-5 20])
% xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
xlabel('\Delta Lapse (LED - Ctrl)')
title(['R: ' num2str(R(2,1)) '  P: ' num2str(P(2,1))])
axis square

subplot(1,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region))& Lapse_sig;
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(D_lapse(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
hold on
idx = strcmp(region,all_region(i_region))& ~Lapse_sig;
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(D_lapse(idx),Temp,'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',[1 1 1],'SizeData',40)

end 
[R,P] = corrcoef(diff(FA,[],2), diff(Lapse,[],2) ); 


ylim([-0.1 0.1])
% xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
xlabel('\Delta Lapse (LED - Ctrl)')
title(['R: ' num2str(R(2,1)) '  P: ' num2str(P(2,1))])
axis square

%% plot delta change over light power seperate by methods

figure
subplot(2,2,1)
for i_b = 1:4
idx = strcmp(Background,B_unique(i_b));
temp = Thresh(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on


end 
ylim([-5 20])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)')
xlabel('light power (mW)')
axis square

subplot(2,2,2)
for i_b = 1:4
idx = strcmp(Background,B_unique(i_b));
temp = FA(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on


end 
ylim([-0.1 0.1])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)')
xlabel('light power (mW)')
axis square

subplot(2,2,3)
for i_b = 1:4
idx = strcmp(Background,B_unique(i_b));
temp = Lapse(idx,:);
Temp = diff(temp,[],2);

scatter(power(idx),Temp,'MarkerEdgeColor',back_color(i_b,:),'MarkerFaceColor',back_color(i_b,:),'SizeData',40)

hold on


end 
ylim([-0.06 0.12])
xlim([0 1])
hline(0,'k:')
set(gca,'TickDir','out','YTick',[-0.06:0.06:0.12])
ylabel('\Delta Lapse (LED - Ctrl)')
xlabel('light power (mW)')
axis square
%% plot percentage change for everything

Colors(1,:) = [0 0 0];
Colors(2:4,:)= lines(3);
all_region = {'V1','LM','AL','PM'};
figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Thresh(idx,:);

exp = sum(idx);
p_change(i_region,1,1) = mean(diff(temp,1,2)./temp(:,1));
p_change(i_region,2,1) = std(diff(temp,1,2)./temp(:,1))./sqrt(exp); 
scatter(repmat(i_region,exp,1),diff(temp,1,2)./temp(:,1),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)./temp(:,1)),std(diff(temp,1,2)./temp(:,1))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')

hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
 text(i_region,19,['n=' num2str(exp)])


end 
ylim([-0.5 1.5])
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta Thresh (LED - Ctrl)/Ctrl')
axis square

subplot(2,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA(idx,:);
exp = sum(idx);
p_change(i_region,1,2) = mean(diff(temp,1,2)./temp(:,1));
p_change(i_region,2,2) = std(diff(temp,1,2)./temp(:,1))./sqrt(exp); 
scatter(repmat(i_region,exp,1),diff(temp,1,2)./temp(:,1),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)./temp(:,1)),std(diff(temp,1,2)./temp(:,1))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
% scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.07,['n=' num2str(exp)])


end 
ylim([-0.8 0.4])
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta FA (LED - Ctrl)/Ctrl')
axis square

subplot(2,2,3)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = d_22(idx,:);
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2)./temp(:,1),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)./temp(:,1)),std(diff(temp,1,2)./temp(:,1))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
% text(i_region,0.7,['n=' num2str(exp)])


end 
ylim([-0.5 0.2])
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta d prime-22.5 (LED - Ctrl)/Ctrl')
axis square

subplot(2,2,4)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = c_22(idx,:);
exp = sum(idx);
scatter(repmat(i_region,exp,1),diff(temp,1,2)./temp(:,1),'MarkerEdgeColor',Colors(i_region,:),'SizeData',10)
hold on
errorbar(i_region,mean(diff(temp,1,2)./temp(:,1)),std(diff(temp,1,2)./temp(:,1))./sqrt(exp),'Color',Colors(i_region,:),'Marker','o')
hold on
%scatter(i_region,mean(diff(temp,1,2)),'MarkerEdgeColor',Colors(i_region,:),'MarkerFaceColor',Colors(i_region,:),'SizeData',40)
%text(i_region,0.7,['n=' num2str(exp)])


end 
ylim([-1 5])
xlim([0 5])
hline(0,'k:')
set(gca,'XTick',1:1:4,'XTickLabel',all_region,'TickDir','out')
ylabel('\Delta c-22.5 (LED - Ctrl)/Ctrl')
axis square

%% save individual plot 
ledcolor = {[0 0 0] [0.2 0.6 1] [1 0 0.6] [0 0.6 0.2]};
for i = [9:12]%1:length(IDs)
    figure
    suptitle(sprintf('%d-%s-%3.2fmW-%s-%s',IDs(i),region{i},power(i),Background{i}))  
    a=5; 
    subplot(1,2,1)
    for i_led = 1:2
        
        Orien = Oriens(i,:,i_led);
        Orien = Orien(~isnan(Orien)); 
        
        maxI = max(Orien);
        minI = min(Orien);
        xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
        h=line(xgrid, fitSall{i,i_led}.modelFun(fitSall{i,i_led}.coefEsts, xgrid), 'Color',ledcolor{i_led});
        hold on;
        plot(fitSall{i,i_led}.intensityX,fitSall{i,i_led}.fractCorrY, 'o','Color',ledcolor{i_led});
        %thresh = coefEsts(1)*[1 1];
        plot(fitSall{i,i_led}.thresh*[1 1], [0 fitSall{i,i_led}.threshY], '--','Color',ledcolor{i_led});
        plot(bootStats{i,i_led}.ci95, fitSall{i,i_led}.threshY*[1 1], 'Color',ledcolor{i_led});
       
        % set limits correctly
        xLim = [min(xgrid) max(xgrid)].* [0.75 1.25];
        xLim = 10.^ceil(log10(xLim) - [1 0]);    
        
    end
    axis([a 100 0 1]) 
    set(gca,'xscale','log','XTick', [10:10:100], 'XTickLabel', {'10';'';'';'';' ';'';' ';' ';'';'100'},'TickDir','out')
    xlabel('Orientation change degree')
    ylabel('Hit Rate')
    axis square
    
    subplot(1,2,2)
    for i_led = 1:2
        
        scatter(0, FA(i,i_led),'MarkerEdgeColor',ledcolor{i_led})
        hold on
        plot([0 0], squeeze(FA_CI(i,i_led,:)), 'Color',ledcolor{i_led})
        
    end
    ylim([0 0.2])
    set(gca,'XTick', [0], 'XTickLabel', {'Collapsed FA'},'TickDir','out')
    ylabel('FA rate')
    axis square
    
    print(sprintf('%d-%s-%3.2fmW-%s.pdf',IDs(i),region{i},power(i),Background{i}),'-dpdf','-fillpage')
    %print(sprintf('%d-%s-%3.2fmW-%s-%s.pdf',IDs(i),region{i},power(i),background{i},virus{i}),'-dpdf','-fillpage')
    close all
end
%% summary of lapse rate, slope,
figure
subplot(1,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = Lapse(idx,:);
exp = sum(idx);
[h_l(i_region),p_l(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region,:),'Marker','o','LineStyle','none')

end 

ylim([0 0.5])
xlim([0.2 4.8])
set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',all_region)%{'V1', 'LM', 'AL','PM'}
ylabel('lapse rate')
axis square

subplot(1,2,2)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = slope(idx,:);
exp = sum(idx);
[h_s(i_region),p_s(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region,:),'Marker','o','LineStyle','none')

end 

ylim([0 3])
xlim([0.2 4.8])
set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',all_region)%{'V1', 'LM', 'AL','PM'}
ylabel('Slope of fitted function')
axis square
%%  summary of reaction time for 90deg, reaction time near threshold orien for hits, reaction time for FAs 


figure
subplot(2,2,1)
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = FA_RT(idx,:)./1000;
exp = sum(idx);
[h_frt(i_region),p_frt(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region,:),'Marker','o','LineStyle','none')

end 

ylim([0.2 0.5])
xlim([0.2 4.8])
set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',all_region)
ylabel('FA-RT (s)')
axis square

subplot(2,2,2)
% RT for easiest target
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = RT_easy(idx,:)./1000;
exp = sum(idx);
[h_ert(i_region),p_ert(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region,:),'Marker','o','LineStyle','none')

end 

ylim([0.2 0.5])
xlim([0.2 4.8])
set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',all_region)
ylabel('RT-easiest target(s)')
axis square
subplot(2,2,3)
% RT for easiest target
for i_region = 1:4
idx = strcmp(region,all_region(i_region));
temp = RT_thresh(idx,:)./1000;
exp = sum(idx);
[h_trt(i_region),p_trt(i_region)]= ttest(temp(:,1),temp(:,2));
plot([i_region-0.3 i_region+0.3],temp, 'color',Colors(i_region,:))

hold on
errorbar([i_region-0.3 i_region+0.3],[mean(temp(:,1)) mean(temp(:,2))],[std(temp(:,1))./sqrt(exp) std(temp(:,2))./sqrt(exp)],'Color',Colors(i_region,:),'Marker','o','LineStyle','none')

end 

ylim([0.2 0.5])
xlim([0.2 4.8])
set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',all_region)
ylabel('RT-threshold ori(s)')
axis square

%% plot threshold and FA rate seperated by preceding trials 
% reorganize it so that Ctrl/Ctrl LED/Ctrl Ctrl/LED LED/LED
Pretrl_FA(:,1) = LEDpre_FA(:,3);
Pretrl_FA(:,2) = LEDpre_FA(:,1);
Pretrl_FA(:,3) = LEDpre_FA(:,4);
Pretrl_FA(:,4) = LEDpre_FA(:,2);

Pretrl_thresh(:,1) = LEDpre_thresh(:,3);
Pretrl_thresh(:,2) = LEDpre_thresh(:,1);
Pretrl_thresh(:,3) = LEDpre_thresh(:,4);
Pretrl_thresh(:,4) = LEDpre_thresh(:,2);
%% thresh
figure
for i_region = 1:4
    subplot(2,2,i_region)
    
    idx = strcmp(region,all_region(i_region));
    temp = Pretrl_thresh(idx,:);
    exp = sum(idx);
    
    plot(1:4,temp, 'color',Colors(i_region,:))
    
    hold on
    errorbar(1:4,mean(temp,1),std(temp,[],1)./sqrt(exp),'Color',Colors(i_region,:),'Marker','o','LineStyle','none')
    text(0.7,48,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 50])
    xlim([0.5 4.5])
    set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',{'Ctrl|Ctrl','LED/Ctrl','Ctrl/LED','LED/LED'})
    ylabel('Threshold (deg)')
    axis square
    xtickangle(30)
end

%% FA
figure
for i_region = 1:4
    subplot(2,2,i_region)
    
    idx = strcmp(region,all_region(i_region));
    temp = Pretrl_FA(idx,:);
    exp = sum(idx);
    
    plot(1:4,temp, 'color',Colors(i_region,:))
    
    hold on
    errorbar(1:4,mean(temp,1),std(temp,[],1)./sqrt(exp),'Color',Colors(i_region,:),'Marker','o','LineStyle','none')
    text(0.7,48,all_region(i_region),'color',Colors(i_region,:))
     ylim([0 0.2])
    xlim([0.5 4.5])
    set(gca,'TickDir','out','XTick',1:1:4, 'XTickLabel',{'Ctrl|Ctrl','LED/Ctrl','Ctrl/LED','LED/LED'})
    ylabel('FA rate')
    axis square
    xtickangle(30)
end
%% stats
i_region =4;
%idx = strcmp(region,all_region(i_region));
%temp = Pretrl_FA(idx,:);
%data = [temp(:,1:2);temp(:,3:4)]; 
idx = strcmp(region,all_region(i_region))& Mice_select';
temp = squeeze(FA_bintrial(idx,:,:)); 
exp = sum(idx);
data = [squeeze(temp(:,:,1))./squeeze(temp(:,1,1)); squeeze(temp(:,:,2))./squeeze(temp(:,1,2))]; 
[p,tbl,stats] = anova2(data,exp);
%c = multcompare(stats,'Estimate','Column'); 
%% compare the p changes across affected areas
P_change = [];
Group = [];
for i_region = 1:3
idx = strcmp(region,all_region(i_region));
temp = c_22(idx,:);
exp = sum(idx);
P_change = [P_change; diff(temp,1,2)./temp(:,1)]; 
Group = [Group; repmat(i_region,exp,1)]; 

end 
[p,tbl,stats] = anova1(P_change,Group); 
%% dependence on trial length both for FA and threshold

figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region));
    temp = squeeze(Thresh_bintrial(idx,:,1));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(Thresh_bintrial(idx,:,2));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(5,48,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 50])
    xlim([1 6])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Thresh')
    axis square
end
figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region));
    temp = squeeze(FA_bintrial(idx,:,1));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(FA_bintrial(idx,:,2));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(1.1,0.19,all_region(i_region),'color',Colors(i_region,:))
    ylim([0 0.2])
    xlim([1 6])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('FA rate')
    axis square
end

%% Only for mice that has enough trials 

figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(Thresh_bintrial(idx,:,1));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(Thresh_bintrial(idx,:,2));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(1.1,48,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    
    ylim([0 50])
    xlim([1 5.5])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Thresh')
    axis square
end

figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(FA_bintrial(idx,:,1));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(FA_bintrial(idx,:,2));
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(1.1,0.19,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    ylim([0 0.2])
    xlim([1 5.5])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('FA rate')
    axis square
end

%% dependence on trial length both for FA and threshold, normalized by the the first bin  

figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region));
    temp = squeeze(Thresh_bintrial(idx,:,1));
    temp = temp./repmat(temp(:,1),1,4);
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(Thresh_bintrial(idx,:,2));
    temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(5,0.9,all_region(i_region),'color',Colors(i_region,:))
    ylim([0.5 1.1])
    xlim([1 6])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Norm. Thresh')
    axis square
end

figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region));
    temp = squeeze(FA_bintrial(idx,:,1));
    temp = temp./repmat(temp(:,1),1,4);
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(FA_bintrial(idx,:,2));
    temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(5,1,all_region(i_region),'color',Colors(i_region,:))
    ylim([0.9 3])
    xlim([1 6])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Norm. FA rate')
    axis square
end

%% for selected mice
figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(Thresh_bintrial(idx,:,1));
    temp = temp./repmat(temp(:,1),1,4);
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(Thresh_bintrial(idx,:,2));
    temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(5,0.9,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    ylim([0.5 1.1])
    xlim([1 5.5])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Norm. Thresh')
    axis square
end

figure
for i_region = 1:4
    
    subplot(2,2,i_region)
    idx = strcmp(region,all_region(i_region))& Mice_select';
    temp = squeeze(FA_bintrial(idx,:,1));
    temp = temp./repmat(temp(:,1),1,4);
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{1}},0)
    hold on
    temp = squeeze(FA_bintrial(idx,:,2));
    temp = temp./repmat(temp(:,1),1,4); % normalized by its own first bin
    shadedErrorBar_ch(cbin./1000, mean(temp,1),std(temp,[],1)./sqrt(sum(idx)),{'Color',ledcolor{2}},0)
    text(5,1,[all_region(i_region) num2str(sum(idx))],'color',Colors(i_region,:))
    ylim([0.9 3])
    xlim([1 5.5])
    set(gca,'TickDir','out')
    xlabel('Trial length (s)')
    ylabel('Norm. FA rate')
    axis square
end
%% stats
data = [];
data = squeeze(FA_bintrial(logical(Mice_select),:,1)); 
%data = data./(data(:,1)); 
[p,tbl,stats] = anova1(data); 
%% get the threshold, FA and lapse value 

uniqueID = unique(IDs); 
for i=1:length(uniqueID)
    idx=[];
    idx = IDs==uniqueID(i); 
    ThreshU(i) = mean(Thresh(idx,1)); 
    FAU(i) = mean(FA(idx,1));
    lapseU(i) = mean(Lapse(idx,1));
    d22U(i)=mean(d_22(idx,1)); 
end 
mean(lapseU)
std(lapseU)./sqrt(length(lapseU))


mean(d22U)
std(d22U)./sqrt(length(d22U))
% mouse = IDs(logical(Mice_select)); 
% data = [];
% data = squeeze(Thresh_bintrial(logical(Mice_select),:,1)); 
% Uni_select = logical([1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 0 0]); 
% %data = data./(data(:,1));
% data = data(Uni_select,:)./(data(Uni_select,1)); 
% [p,tbl,stats] = anova1(data); 

%% chose the unique mice for the reaction time difference
Uni_select = logical([1 1 1 1 1 0 0 1 0 0 1 1 1 1 1 1 1 1 1 1 0]); 
[h,p]=ttest(RT_easy(Uni_select,:), RT_thresh(Uni_select,:)); 


%% register the trial number and sessions
Trial = size(inputcombined.reactTimesMs,2)+ size(inputcombined2.reactTimesMs,2) 
Session = sum(diff(inputcombined.tseq)<0)+1 
%% get reaction time distribution for an example mouse
load('Z:\All_Staff\home\miaomiao\Analysis\Behavior\select and analysis\ISIs\Full duration-led\i527_ISI-LED.mat')
figure
cdfplot([Output{1}.RT.T_RT{1} Output{1}.RT.T_RT{2} Output{1}.RT.T_RT{3}])
hold on
vline([200 550],'k:')
grid off
axis square
xlabel('React time from target onset (ms)')
ylabel('Fraction of trials')
title('example mouse')

