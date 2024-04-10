supp_by_mouse = nan(3,nSess);
facil_by_mouse = nan(3,nSess);

for iSess = 1:nSess
    mouseIndsTemp = mouseInds{iSess};
    red_ind_temp = intersect(red_ind_concat, mouseIndsTemp);
    normDiff_temp = norm_diff(:,:,red_ind_temp);
    N=length(red_ind_temp);

    supp_temp=normDiff_temp(:,:,:)<=-1;
    facil_temp=normDiff_temp(:,:,:)>=1;

    supp_by_mouse(:,iSess)=sum(supp_temp(1,:,:),3)/N;
    facil_by_mouse(:,iSess)=sum(facil_temp(1,:,:),3)/N;

clear mouseIndsTemp red_ind_temp normDiff_temp N supp_temp facil_temp

end

%% summarize across mice
weights=(cellCountsRed(:,1)./length(red_ind_concat))*100;

supp_means= sum((supp_by_mouse.*weights'),2)/100;
supp_std = std(supp_by_mouse,weights',2);
supp_se = supp_std/nSess;

facil_means= sum((facil_by_mouse.*weights'),2)/100;
facil_std = std(facil_by_mouse,weights',2);
facil_se = facil_std/nSess;

%% plotting with weighted means

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_means'],'FaceColor',"#00AFEF",'EdgeColor', [1 1 1]);
hold on
errorbar([1,2,3],supp_means,supp_se,'.k')
xticklabels({'25','50','100'})
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
ylim([0 .4])
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_means'],'FaceColor',"#A8518A",'EdgeColor', [1 1 1]);
hold on
errorbar([1,2,3],facil_means,supp_se,'.k')
xticklabels({'25','50','100'})
title('Facilitated')
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_2C_byMouse_weighted.pdf'),'-dpdf')

%% anova with weighted means

fractSupp = reshape(supp_by_mouse,3*nSess,1);
fractFacil = reshape(facil_by_mouse,3*nSess,1);
contrast = repmat(cons,1,nSess)';
mouseID_col = repelem(unique(mouseID),3,1);
weightsCol = repelem(weights,3,1);

testTable=table(mouseID_col,fractSupp,fractFacil,contrast,weightsCol);

lme_supp_by_corr= fitlme(testTable,'fractSupp~contrast+(1|mouseID_col)', Weights=weightsCol);
anova(lme_supp_by_corr)

lme_facil_by_corr= fitlme(testTable,'fractFacil~contrast+(1|mouseID_col)', Weights=weightsCol);
anova(lme_facil_by_corr)

%% for weakly and strongly correlated cells

supp_by_mouse_low = nan(3,nSess);
facil_by_mouse_low = nan(3,nSess);

for iSess = 1:nSess
    mouseIndsTemp = mouseInds{iSess};
    red_ind_temp = intersect(redLow, mouseIndsTemp);
    normDiff_temp = norm_diff(:,:,red_ind_temp);
    N=length(red_ind_temp);

    supp_temp=normDiff_temp(:,:,:)<=-1;
    facil_temp=normDiff_temp(:,:,:)>=1;

    supp_by_mouse_low(:,iSess)=sum(supp_temp(1,:,:),3)/N;
    facil_by_mouse_low(:,iSess)=sum(facil_temp(1,:,:),3)/N;

clear mouseIndsTemp red_ind_temp normDiff_temp N supp_temp facil_temp

end


supp_by_mouse_high = nan(3,nSess);
facil_by_mouse_high = nan(3,nSess);

for iSess = 1:nSess
    mouseIndsTemp = mouseInds{iSess};
    red_ind_temp = intersect(redHigh, mouseIndsTemp);
    normDiff_temp = norm_diff(:,:,red_ind_temp);
    N=length(red_ind_temp);

    supp_temp=normDiff_temp(:,:,:)<=-1;
    facil_temp=normDiff_temp(:,:,:)>=1;

    supp_by_mouse_high(:,iSess)=sum(supp_temp(1,:,:),3)/N;
    facil_by_mouse_high(:,iSess)=sum(facil_temp(1,:,:),3)/N;

clear mouseIndsTemp red_ind_temp normDiff_temp N supp_temp facil_temp

end
%% plotting with raw means

supp_means_low= mean(supp_by_mouse_low,2);
supp_std_low = std(supp_by_mouse_low,[],2);
supp_se_low = supp_std_low/nSess;

facil_means_low= mean(facil_by_mouse_low,2);
facil_std_low = std(facil_by_mouse_low,[],2);
facil_se_low = facil_std_low/nSess;

supp_means_high= mean(supp_by_mouse_high,2);
supp_std_high = std(supp_by_mouse_high,[],2);
supp_se_high = supp_std_high/nSess;

facil_means_high= mean(facil_by_mouse_high,2);
facil_std_high = std(facil_by_mouse_high,[],2);
facil_se_high = facil_std_high/nSess;

supp_means_low= mean(supp_by_mouse_low,2);
supp_std_low = std(supp_by_mouse_low,[],2);
supp_se_low = supp_std_low/nSess;

facil_means_low= mean(facil_by_mouse_low,2);
facil_std_low = std(facil_by_mouse_low,[],2);
facil_se_low = facil_std_low/nSess;


supp_means_high= mean(supp_by_mouse_high,2);
supp_std_high = std(supp_by_mouse_high,[],2);
supp_se_high = supp_std_high/nSess;

facil_means_high= mean(facil_by_mouse_high,2);
facil_std_high = std(facil_by_mouse_high,[],2);
facil_se_high = facil_std_high/nSess;


figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_means_low'; supp_means_high'],'EdgeColor', [1 1 1]);
hold on
errorbar([.85, 1.85,2.85],supp_means_low,supp_se_low,'.k')
errorbar([1.15, 2.15,3.15],supp_means_high,supp_se_low,'.k')
b(1).FaceColor="#70D0F6";
b(2).FaceColor="#0C8ABB";
xticklabels({'25','50','100'})
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_means_low'; facil_means_high'],'EdgeColor', [1 1 1]);
hold on
errorbar([.85, 1.85,2.85],facil_means_low,supp_se_low,'.k')
errorbar([1.15, 2.15,3.15],facil_means_high,supp_se_low,'.k')
b(1).FaceColor="#C983B1";
b(2).FaceColor="#883367";
xticklabels({'25','50','100'})
title('Facilitated')
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])



%% anova with weighted means
fractSupp = vertcat(reshape(supp_by_mouse_low,3*nSess,1),reshape(supp_by_mouse_high,3*nSess,1));
fractFacil = vertcat(reshape(facil_by_mouse_low,3*nSess,1),reshape(facil_by_mouse_high,3*nSess,1));
contrast = repmat(cons,1,nSess*2)';
corr = repelem(["low","high"],nSess*3)';
mouseID_col = repmat(repelem(unique(mouseID),3,1),2,1);
weightsTemp=RbyExp./sum(RbyExp,2)*100;
weightsCol=repelem(cat(2,weightsTemp(1,:),weightsTemp(2,:))',3);

testTable=table(mouseID_col,fractSupp,fractFacil,corr,contrast,weightsCol);

lme_supp_by_corr= fitlme(testTable,'fractSupp~contrast*corr+(1|mouseID_col)', Weights=weightsCol);
anova(lme_supp_by_corr)

lme_facil_by_corr= fitlme(testTable,'fractFacil~contrast*corr+(1|mouseID_col)', Weights=weightsCol);
anova(lme_facil_by_corr)

%% plotting with weighted means

supp_means_low= sum((supp_by_mouse_low.*weightsTemp(1,:)),2)/100;
supp_std_low = std(supp_by_mouse_low,weightsTemp(1,:),2);
supp_se_low = supp_std_low/nSess;

facil_means_low= sum((facil_by_mouse_low.*weightsTemp(1,:)),2)/100;
facil_std_low = std(facil_by_mouse_low,weightsTemp(1,:),2);
facil_se_low = facil_std_low/nSess;


supp_means_high= sum((supp_by_mouse_high.*weightsTemp(2,:)),2)/100;
supp_std_high = std(supp_by_mouse_high,weightsTemp(2,:),2);
supp_se_high = supp_std_high/nSess;

facil_means_high= sum((facil_by_mouse_high.*weightsTemp(2,:)),2)/100;
facil_std_high = std(facil_by_mouse_high,weightsTemp(2,:),2);
facil_se_high = facil_std_high/nSess;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_means_low'; supp_means_high'],'EdgeColor', [1 1 1]);
hold on
errorbar([.85, 1.85,2.85],supp_means_low,supp_se_low,'.k')
errorbar([1.15, 2.15,3.15],supp_means_high,supp_se_low,'.k')
b(1).FaceColor="#70D0F6";
b(2).FaceColor="#0C8ABB";
xticklabels({'25','50','100'})
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_means_low'; facil_means_high'],'EdgeColor', [1 1 1]);
hold on
errorbar([.85, 1.85,2.85],facil_means_low,supp_se_low,'.k')
errorbar([1.15, 2.15,3.15],facil_means_high,supp_se_low,'.k')
b(1).FaceColor="#C983B1";
b(2).FaceColor="#883367";
xticklabels({'25','50','100'})
title('Facilitated')
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_3C_byMouse_weighted.pdf'),'-dpdf')

%% pairwise ttest, no weighting
[h,p] = ttest(supp_by_mouse_low(3,:),supp_by_mouse_high(3,:))

%% for running vs stationary
supp_by_mouse_stat = nan(3,nSess);
facil_by_mouse_stat = nan(3,nSess);
supp_by_mouse_loc = nan(3,nSess);
facil_by_mouse_loc = nan(3,nSess);

counts = [];

for iSess = 1:nSess
    mouseIndsTemp = mouseInds{iSess};
    red_ind_temp = intersect(red_all, mouseIndsTemp);
    normDiff_temp = norm_diff(:,:,red_ind_temp);
    N=length(red_ind_temp);
    counts=[counts, N];
    supp_temp=normDiff_temp(:,:,:)<=-1;
    facil_temp=normDiff_temp(:,:,:)>=1;

    supp_by_mouse_stat(:,iSess)=sum(supp_temp(1,:,:),3)/N;
    facil_by_mouse_stat(:,iSess)=sum(facil_temp(1,:,:),3)/N;

    supp_by_mouse_loc(:,iSess)=sum(supp_temp(2,:,:),3)/N;
    facil_by_mouse_loc(:,iSess)=sum(facil_temp(2,:,:),3)/N;

clear mouseIndsTemp red_ind_temp normDiff_temp N supp_temp facil_temp

end

%% weighted suppression and facilitation tables
weightedSupp_Stat = (supp_by_mouse_stat .* counts)/length(red_all)

%% get weighted means and se for plotting
weights=counts/length(red_all)*100;
nSess_Temp=sum(counts>0);

supp_means_stat= sum((supp_by_mouse_stat.*weights),2,'omitmissing')/100;
supp_std_stat = std(supp_by_mouse_stat,weights,2,'omitmissing');
supp_se_stat = supp_std_stat/nSess_Temp;

facil_means_stat= sum((facil_by_mouse_stat.*weights),2,'omitmissing')/100;
facil_std_stat = std(facil_by_mouse_stat,weights,2,'omitmissing');
facil_se_stat = facil_std_stat/nSess_Temp;


supp_means_loc= sum((supp_by_mouse_loc.*weights),2,'omitmissing')/100;
supp_std_loc = std(supp_by_mouse_loc,weights,2,'omitmissing');
supp_se_loc = supp_std_loc/nSess_Temp;

facil_means_loc= sum((facil_by_mouse_loc.*weights),2,'omitmissing')/100;
facil_std_loc = std(facil_by_mouse_loc,weights,2,'omitmissing');
facil_se_loc = facil_std_loc/nSess_Temp;

%% plotting
figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_means_stat'; supp_means_loc'],'EdgeColor', [1 1 1]);
hold on
errorbar([.85, 1.85,2.85],supp_means_stat,supp_se_stat,'.k')
errorbar([1.15, 2.15,3.15],supp_means_loc,supp_se_stat,'.k')
b(1).FaceColor="#70D0F6";
b(2).FaceColor="#0C8ABB";
xticklabels({'25','50','100'})
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_means_stat'; facil_means_loc'],'EdgeColor', [1 1 1]);
hold on
errorbar([.85, 1.85,2.85],facil_means_stat,supp_se_stat,'.k')
errorbar([1.15, 2.15,3.15],facil_means_loc,supp_se_stat,'.k')
b(1).FaceColor="#C983B1";
b(2).FaceColor="#883367";
xticklabels({'25','50','100'})
title('Facilitated')
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_4C_byMouse_weighted.pdf'),'-dpdf')
%% anova with weighted means
fractSupp = vertcat(reshape(supp_by_mouse_stat,3*nSess,1),reshape(supp_by_mouse_loc,3*nSess,1));
fractFacil = vertcat(reshape(facil_by_mouse_stat,3*nSess,1),reshape(facil_by_mouse_loc,3*nSess,1));
contrast = repmat(cons,1,nSess*2)';
state = repelem(["stat","loc"],nSess*3)';
mouseID_col = repmat(repelem(unique(mouseID),3,1),2,1);
weightsCol= repmat(repelem(weights,1,3)',2,1);
testTable=table(mouseID_col,fractSupp,fractFacil,contrast,state,weightsCol);

lme_supp_by_corr= fitlme(testTable,'fractSupp~contrast*state+(1|mouseID_col)', Weights=weightsCol);
anova(lme_supp_by_corr)

lme_facil_by_corr= fitlme(testTable,'fractFacil~contrast*state+(1|mouseID_col)', Weights=weightsCol);
anova(lme_facil_by_corr)


%% pairwise ttest, no weighting
[h1,p1] = ttest(facil_by_mouse_stat(1,:),facil_by_mouse_loc(1,:));
[h2,p2] = ttest(facil_by_mouse_stat(2,:),facil_by_mouse_loc(2,:));

[p1*2, p2*2]

[h1,p1] = ttest(supp_by_mouse_stat(1,:),supp_by_mouse_loc(1,:));
[h2,p2] = ttest(supp_by_mouse_stat(3,:),supp_by_mouse_loc(2,:));

[p1*2, p2*2]