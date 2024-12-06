%import avg response to each ori for each cell on each day

a = d1_fits.avgResponseEaOri(match_all,:);
b = d2_fits.avgResponseEaOri(match_all,:);
c = d3_fits.avgResponseEaOri(match_all,:);
d = d4_fits.avgResponseEaOri(match_all,:);

%%
%make correlation matrices of base1 to other sessions for avg ori
basecorr = corr(a,b);
figure;
subplot(2,2,1)
heatmap(basecorr)
xlabel('Orientation Base1')
ylabel('Orientation Base2')
colormap default
subplot(2,2,2)
darkcorr = corr(a,c);
heatmap(darkcorr)
xlabel('Orientation Base1')
ylabel('Orientation Post-dark')
colormap default
subplot(2,2,3)
lightcorr = corr(a,d);
heatmap(lightcorr)
xlabel('Orientation Base1')
ylabel('Orientation Post-dark+7d')
colormap default

%%
%make and plot diagonals (same ori across sessions)
new_mat = [];
new_mat(:,1) = diag(basecorr);
new_mat(:,2) = diag(darkcorr);
new_mat(:,3) = diag(lightcorr);

figure;
plot(new_mat)
xlim([0,9])
ylim([0,1])
xlabel('Orientation')
ylabel('Population Vector Correlation')
legend(['Base1 - Base2'], ['Base1 - Post-dark'], ['Base1 - Post-dark+7d'], 'Location', 'Best')
newcolors = [0 0 0
             .5 .5 .5
             .8 .8 .8];
colororder(newcolors)

%%
%plot mean across all oris between sessions
figure;
errorbar(mean(new_mat), (std(new_mat)/sqrt(length(new_mat))))
xlim([0 4])
ylim([0 1])

%%
%base1 correlated with itself
base1corr_base1corr = corr(a,a);
off_diag_base1_base1 = [];

for i = 1:(length(base1corr_base1corr)-1)
    off_diag_base1_base1(i) = base1corr_base1corr(i,i+1);
end
off_diag_base1_base1 = off_diag_base1_base1';

%base2 correlated w/ itself
base2corr_base2corr = corr(b,b);
off_diag_base2_base2 = [];

for i = 1:(length(base2corr_base2corr)-1)
    off_diag_base2_base2(i) = base2corr_base2corr(i,i+1);
end
off_diag_base2_base2 = off_diag_base2_base2';

%post-dark correlated w/ itself
darkcorr_darkcorr = corr(c,c);
off_diag_dark_dark = [];

for i = 1:(length(darkcorr_darkcorr)-1)
    off_diag_dark_dark(i) = darkcorr_darkcorr(i,i+1);
end
off_diag_dark_dark = off_diag_dark_dark';

%postdark +7d light correlated w/ itself
lightcorr_lightcorr = corr(d,d);
off_diag_light_light = [];

for i = 1:(length(lightcorr_lightcorr)-1)
    off_diag_light_light(i) = lightcorr_lightcorr(i,i+1);
end
off_diag_light_light = off_diag_light_light';


%%
new_mat_off_diag = [];
new_mat_off_diag(:,1) = off_diag_base1_base1;
new_mat_off_diag(:,2) = off_diag_base2_base2;
new_mat_off_diag(:,3) = off_diag_dark_dark;
new_mat_off_diag(:,4) = off_diag_light_light;


figure;
errorbar(mean(new_mat_off_diag), (std(new_mat_off_diag)/sqrt(length(new_mat_off_diag))))
xlim([0 5])
ylim([0 1])


%%
figure;
errorbar(mean(new_mat), (std(new_mat)/sqrt(length(new_mat))))
hold on;
errorbar(mean(new_mat_off_diag), (std(new_mat_off_diag)/sqrt(length(new_mat_off_diag))))


figure;
errorbar([2,4,6], mean(new_mat), (std(new_mat)/sqrt(length(new_mat))), '-o')
hold on;
errorbar([1,3,5,7], mean(new_mat_off_diag), (std(new_mat_off_diag)/sqrt(length(new_mat_off_diag))), '-o')
xlim([0 8])
xticks([0:8])
xticklabels({'','Base1 w/self','Base1 to Base2','Base2 w/self','Base1 to Post-dark','Post-dark w/self','Base1 to Post-dark+7d','Post-dark+7d w/self',''})
legend(['Between-session (Same ori)'], ['Within-session (n-1 ori)'], 'Location', 'Best')
ylabel('Population Correlation')
print(fullfile(newfnout, ['test', '_test.pdf']), '-dpdf', '-bestfit')


%%

%%
%doing the same as above but each session is compared to the previous
%(rather than to session 1)

%%
figure;
subplot(2,2,1)
heatmap(corr(a,b))
xlabel('Orientation Base1')
ylabel('Orientation Base2')
colormap default
subplot(2,2,2)
heatmap(corr(b,c))
xlabel('Orientation Base2')
ylabel('Orientation Post-dark')
colormap default
subplot(2,2,3)
heatmap(corr(c,d))
xlabel('Orientation Post-dark')
ylabel('Orientation Post-dark+7d')
colormap default




%%

%make and plot diagonals (same ori across sessions)
new_mat_2 = [];
new_mat_2(:,1) = diag(corr(a,b));
new_mat_2(:,2) = diag(corr(b,c));
new_mat_2(:,3) = diag(corr(c,d));

%%
figure;
errorbar([2,4,6], mean(new_mat_2), (std(new_mat_2)/sqrt(length(new_mat_2))), '-o')
hold on;
errorbar([1,3,5,7], mean(new_mat_off_diag), (std(new_mat_off_diag)/sqrt(length(new_mat_off_diag))), '-o')
xlim([0 8])
xticks([0:8])
xticklabels({'','Base1 w/self','Base1 to Base2','Base2 w/self','Base2 to Post-dark','Post-dark w/self','Post-dark to Post-dark+7d','Post-dark+7d w/self',''})
legend(['Between-session (Same ori)'], ['Within-session (n-1 ori)'], 'Location', 'Best')
ylabel('Population Correlation')


%%
%farther off-diagonals
%real one!!!


%base1_base1
off_diag_base1_base1_22point5diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 8
        off_diag_base1_base1_22point5diff(i) = base1corr_base1corr(i-7,i);
    else
        off_diag_base1_base1_22point5diff(i) = base1corr_base1corr(i+1,i);
    end
       
end
off_diag_base1_base1_22point5diff = off_diag_base1_base1_22point5diff';


off_diag_base1_base1_45diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 7
        off_diag_base1_base1_45diff(i) = base1corr_base1corr(i-6,i);
    elseif i == 8
        off_diag_base1_base1_45diff(i) = base1corr_base1corr(i-6,i);
    else
        off_diag_base1_base1_45diff(i) = base1corr_base1corr(i+2,i);
    end
       
end
off_diag_base1_base1_45diff = off_diag_base1_base1_45diff';


off_diag_base1_base1_67point5diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 6
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i-5,i);
    elseif i == 7
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i-5,i);
    elseif i == 8
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i-5,i);
    else
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i+3,i);
    end
       
end
off_diag_base1_base1_67point5diff = off_diag_base1_base1_67point5diff';


off_diag_base1_base1_90diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 5
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    elseif i == 6
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    elseif i == 7
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    elseif i == 8
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    else
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i+4,i);
    end
       
end
off_diag_base1_base1_90diff = off_diag_base1_base1_90diff';

%base2_base2
off_diag_base2_base2_22point5diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 8
        off_diag_base2_base2_22point5diff(i) = base2corr_base2corr(i-7,i);
    else
        off_diag_base2_base2_22point5diff(i) = base2corr_base2corr(i+1,i);
    end
       
end
off_diag_base2_base2_22point5diff = off_diag_base2_base2_22point5diff';


off_diag_base2_base2_45diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 7
        off_diag_base2_base2_45diff(i) = base2corr_base2corr(i-6,i);
    elseif i == 8
        off_diag_base2_base2_45diff(i) = base2corr_base2corr(i-6,i);
    else
        off_diag_base2_base2_45diff(i) = base2corr_base2corr(i+2,i);
    end
       
end
off_diag_base2_base2_45diff = off_diag_base2_base2_45diff';


off_diag_base2_base2_67point5diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 6
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i-5,i);
    elseif i == 7
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i-5,i);
    elseif i == 8
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i-5,i);
    else
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i+3,i);
    end
       
end
off_diag_base2_base2_67point5diff = off_diag_base2_base2_67point5diff';


off_diag_base2_base2_90diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 5
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    elseif i == 6
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    elseif i == 7
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    elseif i == 8
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    else
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i+4,i);
    end
       
end
off_diag_base2_base2_90diff = off_diag_base2_base2_90diff';

%dark_dark
off_diag_dark_dark_22point5diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 8
        off_diag_dark_dark_22point5diff(i) = darkcorr_darkcorr(i-7,i);
    else
        off_diag_dark_dark_22point5diff(i) = darkcorr_darkcorr(i+1,i);
    end
       
end
off_diag_dark_dark_22point5diff = off_diag_dark_dark_22point5diff';


off_diag_dark_dark_45diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 7
        off_diag_dark_dark_45diff(i) = darkcorr_darkcorr(i-6,i);
    elseif i == 8
        off_diag_dark_dark_45diff(i) = darkcorr_darkcorr(i-6,i);
    else
        off_diag_dark_dark_45diff(i) = darkcorr_darkcorr(i+2,i);
    end
       
end
off_diag_dark_dark_45diff = off_diag_dark_dark_45diff';


off_diag_dark_dark_67point5diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 6
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i-5,i);
    elseif i == 7
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i-5,i);
    elseif i == 8
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i-5,i);
    else
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i+3,i);
    end
       
end
off_diag_dark_dark_67point5diff = off_diag_dark_dark_67point5diff';


off_diag_dark_dark_90diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 5
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    elseif i == 6
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    elseif i == 7
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    elseif i == 8
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    else
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i+4,i);
    end
       
end
off_diag_dark_dark_90diff = off_diag_dark_dark_90diff';

%light_light
off_diag_light_light_22point5diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 8
        off_diag_light_light_22point5diff(i) = lightcorr_lightcorr(i-7,i);
    else
        off_diag_light_light_22point5diff(i) = lightcorr_lightcorr(i+1,i);
    end
       
end
off_diag_light_light_22point5diff = off_diag_light_light_22point5diff';


off_diag_light_light_45diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 7
        off_diag_light_light_45diff(i) = lightcorr_lightcorr(i-6,i);
    elseif i == 8
        off_diag_light_light_45diff(i) = lightcorr_lightcorr(i-6,i);
    else
        off_diag_light_light_45diff(i) = lightcorr_lightcorr(i+2,i);
    end
       
end
off_diag_light_light_45diff = off_diag_light_light_45diff';


off_diag_light_light_67point5diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 6
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i-5,i);
    elseif i == 7
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i-5,i);
    elseif i == 8
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i-5,i);
    else
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i+3,i);
    end
       
end
off_diag_light_light_67point5diff = off_diag_light_light_67point5diff';


off_diag_light_light_90diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 5
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    elseif i == 6
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    elseif i == 7
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    elseif i == 8
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    else
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i+4,i);
    end
       
end
off_diag_light_light_90diff = off_diag_light_light_90diff';


%%

new_mat_off_diag_22point5diff(:,1) = off_diag_base1_base1_22point5diff;
new_mat_off_diag_22point5diff(:,2) = off_diag_base2_base2_22point5diff;
new_mat_off_diag_22point5diff(:,3) = off_diag_dark_dark_22point5diff;
new_mat_off_diag_22point5diff(:,4) = off_diag_light_light_22point5diff;

new_mat_off_diag_45diff(:,1) = off_diag_base1_base1_45diff;
new_mat_off_diag_45diff(:,2) = off_diag_base2_base2_45diff;
new_mat_off_diag_45diff(:,3) = off_diag_dark_dark_45diff;
new_mat_off_diag_45diff(:,4) = off_diag_light_light_45diff;

new_mat_off_diag_67point5diff(:,1) = off_diag_base1_base1_67point5diff;
new_mat_off_diag_67point5diff(:,2) = off_diag_base2_base2_67point5diff;
new_mat_off_diag_67point5diff(:,3) = off_diag_dark_dark_67point5diff;
new_mat_off_diag_67point5diff(:,4) = off_diag_light_light_67point5diff;

new_mat_off_diag_90diff(:,1) = off_diag_base1_base1_90diff;
new_mat_off_diag_90diff(:,2) = off_diag_base2_base2_90diff;
new_mat_off_diag_90diff(:,3) = off_diag_dark_dark_90diff;
new_mat_off_diag_90diff(:,4) = off_diag_light_light_90diff;

%%
figure;
errorbar([2,4,6,8], mean(new_mat_off_diag_22point5diff), (std(new_mat_off_diag_22point5diff)/sqrt(length(new_mat_off_diag_22point5diff))), '-o', 'Color', [0 0 0], "MarkerSize", 9)
hold on
errorbar([2,4,6,8], mean(new_mat_off_diag_45diff), (std(new_mat_off_diag_45diff)/sqrt(length(new_mat_off_diag_45diff))), '-diamond', 'Color', [.3 .3 .3], "MarkerSize", 9)
errorbar([2,4,6,8], mean(new_mat_off_diag_67point5diff), (std(new_mat_off_diag_67point5diff)/sqrt(length(new_mat_off_diag_67point5diff))), '-^', 'Color', [.6 .6 .6], "MarkerSize", 9)
errorbar([2,4,6,8], mean(new_mat_off_diag_90diff), (std(new_mat_off_diag_90diff)/sqrt(length(new_mat_off_diag_90diff))), '-square', 'Color', [.9 .9 .9], "MarkerSize", 10)
xlim([0 10])
xticks([0:10])
xticklabels({'','','Base1 w/self','','Base2 w/self','','Post-dark w/self','','Post-dark+7d w/self'})
legend(['22.5 degrees away'], ['45 degrees away'], ['67.5 degrees away'],  ['90 degrees away'], 'FontSize', 12, 'Location', 'Best')
ylabel('Population Correlation')
print(fullfile(newfnout, [mouse,  '_ori_drift_test.pdf']), '-dpdf', '-bestfit')


%%
figure;
errorbar([2,4,6], mean(new_mat_2), (std(new_mat_2)/sqrt(length(new_mat_2))), '-o', 'Color', 'red', "MarkerSize", 9)
hold on;
errorbar([1,3,5,7], mean(new_mat_off_diag_22point5diff), (std(new_mat_off_diag_22point5diff)/sqrt(length(new_mat_off_diag_22point5diff))), '-o', 'Color', [0 0 0], "MarkerSize", 9)
errorbar([1,3,5,7], mean(new_mat_off_diag_45diff), (std(new_mat_off_diag_45diff)/sqrt(length(new_mat_off_diag_45diff))), '-diamond', 'Color', [.3 .3 .3], "MarkerSize", 9)
errorbar([1,3,5,7], mean(new_mat_off_diag_67point5diff), (std(new_mat_off_diag_67point5diff)/sqrt(length(new_mat_off_diag_67point5diff))), '-^', 'Color', [.6 .6 .6], "MarkerSize", 9)
errorbar([1,3,5,7], mean(new_mat_off_diag_90diff), (std(new_mat_off_diag_90diff)/sqrt(length(new_mat_off_diag_90diff))), '-square', 'Color', [.9 .9 .9], "MarkerSize", 10)
xlim([0 8])
xticks([0:8])
xticklabels({'','Base1 w/self','Base1 to Base2','Base2 w/self','Base2 to Post-dark','Post-dark w/self','Post-dark to Post-dark+7d','Post-dark+7d w/self',''})
legend(['Same ori (b/w session)'],['22.5 degrees away (w/i session)'], ['45 degrees away (w/i session)'], ['67.5 degrees away (w/i session)'],  ['90 degrees away (w/i session)'], 'FontSize', 12, 'Location', 'Best')
ylabel('Population Correlation')


%%



%%






%%
%off diagonals across sessions (DO NOT USE BELOW THIS LINE)
off_diag_base_mat = [];

for i = 1:length(basecorr)
    if i == 1
        off_diag_base_mat(i,i) = mean([basecorr(i,i+1), basecorr(i+1,i)])
    elseif i == 8
        off_diag_base_mat(i,i) = mean([basecorr(i,i-1), basecorr(i-1,i)])
    else
        off_diag_base_mat(i,i) = mean([basecorr(i,i-1), basecorr(i-1,i), basecorr(i,i+1), basecorr(i+1,i)])
    end
end


off_diag_dark_mat = [];

for i = 1:length(darkcorr)
    if i == 1
        off_diag_dark_mat(i,i) = mean([darkcorr(i,i+1), darkcorr(i+1,i)])
    elseif i == 8
        off_diag_dark_mat(i,i) = mean([darkcorr(i,i-1), darkcorr(i-1,i)])
    else
        off_diag_dark_mat(i,i) = mean([darkcorr(i,i-1), darkcorr(i-1,i), darkcorr(i,i+1), darkcorr(i+1,i)])
    end
end


off_diag_light_mat = [];

for i = 1:length(lightcorr)
    if i == 1
        off_diag_light_mat(i,i) = mean([lightcorr(i,i+1), lightcorr(i+1,i)])
    elseif i == 8
        off_diag_light_mat(i,i) = mean([lightcorr(i,i-1), lightcorr(i-1,i)])
    else
        off_diag_light_mat(i,i) = mean([lightcorr(i,i-1), lightcorr(i-1,i), lightcorr(i,i+1), lightcorr(i+1,i)])
    end
end

%%
new_mat_off_diag = [];
new_mat_off_diag(:,1) = diag(off_diag_base_mat)
new_mat_off_diag(:,2) = diag(off_diag_dark_mat)
new_mat_off_diag(:,3) = diag(off_diag_light_mat)


figure;
plot(new_mat_off_diag)
xlim([0,9])
ylim([0,1])
xlabel('Orientation')
ylabel('Population Vector Correlation')
legend(['Base1 - Base2'], ['Base1 - Post-dark'], ['Base1 - Post-dark+7d'], 'Location', 'Best')
newcolors = [0 0 0
             .5 .5 .5
             .8 .8 .8];
colororder(newcolors)

figure;
errorbar(mean(new_mat_off_diag), (std(new_mat_off_diag)/sqrt(length(new_mat_off_diag))))
xlim([0 4])
ylim([0 1])

%%

figure;
errorbar(mean(new_mat), (std(new_mat)/sqrt(length(new_mat))))
xlim([0 4])
ylim([0 1])
hold on
errorbar(mean(new_mat_off_diag), (std(new_mat_off_diag)/sqrt(length(new_mat_off_diag))))
xlim([0 4])
ylim([0 1])
xlabel('Sessions')
ylabel('Population Vector Correlation')
legend(['Same Orientation'], ['1-Off Orientation'], 'Location', 'Best')
xticklabels({'', 'Base1-Base2', 'Base1-Dark', 'Base1-Dark+7d', ''})


