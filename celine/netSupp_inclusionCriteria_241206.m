%find cells that are well centered
%centered = find(dists_concat<10);  
centered = find(h_concat(:,1,nCons)==1); %here using the definition that must respond to small high con stim
largeSizeHighCon = find(h_concat(:,nSizes,nCons)==1);
% cells responsive at the smallest size and highest contrast are "centered"
centerIN = intersect(centered,interNrns);
centerPyr = intersect(centered,pyrCells);
%% timecourses with half rise and half decay averages for cells responsive at each condition


[n n2] = subplotn(nSizes*nCons);
x=1;

figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
      
        responsiveTheseTrials = find(h_concat(:,iSize,iCon));
        responsiveCentered = intersect(responsiveTheseTrials, centered);
        thesePyr=intersect(pyrCells,responsiveCentered);
        theseIN=intersect(interNrns,responsiveCentered);

        temp_mean1 = mean(TCs_stat(:,iSize,iCon,thesePyr),4,"omitnan");
        temp_se1 = std(TCs_stat(:,iSize,iCon,thesePyr),[],4,"omitnan")/sqrt(length(thesePyr));

        temp_mean2 = mean(TCs_stat(:,iSize,iCon,theseIN),4,"omitnan");
        temp_se2 = std(TCs_stat(:,iSize,iCon,theseIN),[],4,"omitnan")/sqrt(length(theseIN));

        mean_halfPeak_Pyr=mean(peak_time(iSize,iCon,thesePyr,3),3,"omitmissing");
        mean_halfDecay_Pyr=(mean(peak_time(iSize,iCon,thesePyr,4),3,"omitmissing"));

        mean_halfPeak_IN=mean(peak_time(iSize,iCon,theseIN,3),3,"omitmissing");
        mean_halfDecay_IN=(mean(peak_time(iSize,iCon,theseIN,4),3,"omitmissing"));

        subplot(n,n2,x)
        vline(mean_halfPeak_Pyr,'--b')
        vline(mean_halfDecay_Pyr,'--m')
        vline(mean_halfPeak_IN,'b')
        vline(mean_halfDecay_IN,'m')
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'g');
        hold on
        shadedErrorBar(t,temp_mean1,temp_se1);
        hold on
        alpha(.5)
        hold on
        ylim([-.03 .1])
        xlim([-.1 .3])
        
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
x0=1;
y0=1;
width=7;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Stationary')

sgtitle('For cells responsive at each condition (within 7.5 deg responders)')
print('TCs_respAtEach_7.5deg.pdf', '-dpdf');
%% timecourses with half rise and half decay averages for centered cells


[n n2] = subplotn(nSizes*nCons);
x=1;

figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
     

        temp_mean1 = mean(TCs_stat(:,iSize,iCon,centerPyr),4,"omitnan");
        temp_se1 = std(TCs_stat(:,iSize,iCon,centerPyr),[],4,"omitnan")/sqrt(length(centerPyr));

        temp_mean2 = mean(TCs_stat(:,iSize,iCon,centerIN),4,"omitnan");
        temp_se2 = std(TCs_stat(:,iSize,iCon,centerIN),[],4,"omitnan")/sqrt(length(centerIN));

        mean_halfPeak_Pyr=mean(peak_time(iSize,iCon,centerPyr,3),3,"omitmissing");
        mean_halfDecay_Pyr=(mean(peak_time(iSize,iCon,centerPyr,4),3,"omitmissing"));

        mean_halfPeak_IN=mean(peak_time(iSize,iCon,centerIN,3),3,"omitmissing");
        mean_halfDecay_IN=(mean(peak_time(iSize,iCon,centerIN,4),3,"omitmissing"));

        subplot(n,n2,x)
        vline(mean_halfPeak_Pyr,'--b')
        vline(mean_halfDecay_Pyr,'--m')
        vline(mean_halfPeak_IN,'b')
        vline(mean_halfDecay_IN,'m')
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'g');
        hold on
        shadedErrorBar(t,temp_mean1,temp_se1);
        hold on
        alpha(.5)
        hold on
        ylim([-.03 .1])
        xlim([-.1 .3])
        
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
x0=1;
y0=1;
width=7;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Stationary')

sgtitle('RF dist < 10, responsive to any stim')
%sgtitle('Responsive to small size high con')
print('TCs_centered_1S.pdf', '-dpdf');

%% rise time, decay time, and FWHM for cells responsive at each size, 80% contrast
temp_mean1 = nan(nCons,nSizes,3);
temp_se1 = nan(nCons,nSizes,3);
temp_mean2 = nan(nCons,nSizes,3);
temp_se2 = nan(nCons,nSizes,3);
IN_ns = [];
Pyr_ns = [];
for iCon = 4:nCons
    for iSize = 1:nSizes %loop through the sizes
        
        
        responsiveTheseTrials = find(h_concat(:,iSize,iCon));
        responsiveCentered = intersect(responsiveTheseTrials, centered);
        thesePyr=intersect(pyrCells,responsiveCentered);
        theseIN=intersect(interNrns,responsiveCentered);

        temp_mean1(iCon,iSize,1:2) = mean(peak_time(iSize,iCon,thesePyr,3:4),3,"omitnan");
        temp_se1(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,thesePyr,3:4),[],3,"omitnan"))./length(thesePyr);
        temp_mean1(iCon,iSize,3)=mean(fwhm(iSize,iCon,thesePyr),3,"omitnan");
        temp_se1(iCon,iSize,3)=(std(fwhm(iSize,iCon,thesePyr),[],3,"omitnan"))./length(thesePyr);


        temp_mean2(iCon,iSize,1:2) =  mean(peak_time(iSize,iCon,theseIN,3:4),3,"omitnan");
        temp_se2(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,theseIN,3:4),[],3,"omitnan"))./length(theseIN);
        temp_mean2(iCon,iSize,3)=mean(fwhm(iSize,iCon,theseIN),3,"omitnan");
        temp_se2(iCon,iSize,3)=(std(fwhm(iSize,iCon,theseIN),[],3,"omitnan"))./length(theseIN);

IN_ns = [IN_ns, length(theseIN)];
Pyr_ns = [Pyr_ns, length(thesePyr)];
    end
end

figure;
subplot(1,3,1)
errorbar(Sizes,temp_mean2(4,:,1),temp_se2(4,:,1),'Color',	"#4e701f",'LineStyle','none','Marker','o');
hold on
errorbar(Sizes,temp_mean1(4,:,1),temp_se1(4,:,1),'k','LineStyle','none','Marker','o');
ylim([0,.2])
xticks(Sizes)
xlabel('Size')
ylabel('Seconds')
title('Rise time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,2)
errorbar(Sizes,temp_mean2(4,:,2),temp_se2(4,:,2),'Color',	"#4e701f",'LineStyle','none','Marker','o');
hold on
errorbar(Sizes,temp_mean1(4,:,2),temp_se1(4,:,2),'k','LineStyle','none','Marker','o');
ylim([.1,.3])
xticks(Sizes)
xlabel('Size')
title('Decay time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,3)
errorbar(Sizes,temp_mean2(4,:,3),temp_se2(4,:,3),'Color',	"#4e701f",'LineStyle','none','Marker','o');
hold on
errorbar(Sizes,temp_mean1(4,:,3),temp_se1(4,:,3),'k','LineStyle','none','Marker','o');
ylim([0,.2])
xticks(Sizes)
xlabel('Size')
title('FWHM')
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=7;
height=2.4;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Responsive at each (within 7.5 deg responders)')

print('dynamics_respAtEach_7.5deg', '-dpdf');
%% rise time, decay time, and FWHM for all centered cells 
temp_mean1 = nan(nCons,nSizes,3);
temp_se1 = nan(nCons,nSizes,3);
temp_mean2 = nan(nCons,nSizes,3);
temp_se2 = nan(nCons,nSizes,3);

for iCon = 1:nCons
    for iSize = 1:nSizes %loop through the sizes

        temp_mean1(iCon,iSize,1:2) = mean(peak_time(iSize,iCon,centerPyr,3:4),3,"omitnan");
        temp_se1(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,centerPyr,3:4),[],3,"omitnan"))./length(centerPyr);
        temp_mean1(iCon,iSize,3)=mean(fwhm(iSize,iCon,centerPyr),3,"omitnan");
        temp_se1(iCon,iSize,3)=(std(fwhm(iSize,iCon,centerPyr),[],3,"omitnan"))./length(centerPyr);


        temp_mean2(iCon,iSize,1:2) =  mean(peak_time(iSize,iCon,centerIN,3:4),3,"omitnan");
        temp_se2(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,centerIN,3:4),[],3,"omitnan"))./length(centerIN);
        temp_mean2(iCon,iSize,3)=mean(fwhm(iSize,iCon,centerIN),3,"omitnan");
        temp_se2(iCon,iSize,3)=(std(fwhm(iSize,iCon,centerIN),[],3,"omitnan"))./length(centerIN);


    end
end

figure;
subplot(1,3,1)
errorbar(Sizes,temp_mean2(4,:,1),temp_se2(4,:,1),'Color',	"#4e701f",'Marker','o');
hold on
errorbar(Sizes,temp_mean1(4,:,1),temp_se1(4,:,1),'k','Marker','o');
ylim([0,.2])
xticks(Sizes)
xlabel('Size')
ylabel('Seconds')
title('Rise time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,2)
errorbar(Sizes,temp_mean2(4,:,2),temp_se2(4,:,2),'Color',	"#4e701f",'Marker','o');
hold on
errorbar(Sizes,temp_mean1(4,:,2),temp_se1(4,:,2),'k','Marker','o');
ylim([.1,.3])
xticks(Sizes)
xlabel('Size')
title('Decay time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,3)
errorbar(Sizes,temp_mean2(4,:,3),temp_se2(4,:,3),'Color',	"#4e701f",'Marker','o');
hold on
errorbar(Sizes,temp_mean1(4,:,3),temp_se1(4,:,3),'k','Marker','o');
ylim([0,.2])
xticks(Sizes)
xlabel('Size')
title('FWHM')
set(gca, 'TickDir', 'out')
box off

sgtitle('RF dist < 10, responsive to any stim')

x0=5;
y0=5;
width=7;
height=2.4;
set(gcf,'units','inches','position',[x0,y0,width,height])
