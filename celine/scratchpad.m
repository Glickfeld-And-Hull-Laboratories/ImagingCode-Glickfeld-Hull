
figure;
subplot(1,2,1)
boxchart(squeeze(norm_diff(1,:,redLow))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,redLow))',20,[.79 .25 .32], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)

xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylabel('Normalized difference')
ylim([-5 5])
title('Weakly correlated')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart(squeeze(norm_diff(1,:,redHigh))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,redHigh))',20,[.79 .25 .32], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylim([-5 5])
title('Strongly correlated')
hold off
set(gca,'TickDir','out')
box off
x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% original means
og_means = [mean(squeeze(norm_diff(1,2,redLow)),'omitmissing'), mean(squeeze(norm_diff(1,2,redHigh)),'omitmissing')];
%%

%norm_diff_red = norm_diff(:,:,red_all);

facil_red_low_=norm_diff(:,:,redLow)>=1;
supp_red_low=norm_diff(:,:,redLow)<=-1;
facil_red_high_=norm_diff(:,:,redHigh)>=1;
supp_red_high=norm_diff(:,:,redHigh)<=-1;


N1=length(redLow);
N2=length(redHigh);
facil_table_low = sum(facil_red_low_(1,:,:),3)/N1;
supp_table_low = sum(supp_red_low(1,:,:),3)/N1;
facil_table_high = sum(facil_red_high_(1,:,:),3)/N2;
supp_table_high = sum(supp_red_high(1,:,:),3)/N2;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_low; supp_table_high],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
ylim([0 .6])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_low; facil_table_high],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .6])
title('Facilitated')
ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off
