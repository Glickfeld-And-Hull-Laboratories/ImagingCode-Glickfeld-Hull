%

%[h,p1,ci,stats1] = vartestn(norm_diff(2,1,red_all),norm_diff(2,2,red_all),'Tail','left'); %25 vs 50
%Make sumary table
outputColStat=(squeeze(norm_diff(1,:,red_all)));
outputColStat=reshape(outputColStat,[],1);
outputColLoc=(squeeze(norm_diff(1,:,red_all)));
outputColLoc=reshape(outputColLoc,[],1);
behStateCol = repelem(["stat" "loc"],1,((3*length(red_all))))';
contrastCol=repmat(cons,1,(2*length(red_all)))';
outputCol = vertcat(outputColStat,outputColLoc);

normDiff_summary = table(outputCol,behStateCol,contrastCol, ...
    'VariableNames',{'normDiff' 'behStat' 'contrast'});
clear outputCol behStateCol contrastCol outputColStat outputColLoc

%%
vartestn(normDiff_summary.normDiff, normDiff_summary.contrast,'TestType','Bartlett')
vartestn(normDiff_summary.normDiff, normDiff_summary.behStat,'TestType','Bartlett')
%%



%compute chi squares for facilitation, low vs. high
%25

n1=facil_table_stat_low(1)*N1;
n2=facil_table_stat_high(1)*N2;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=facil_table_stat_low(2)*N1;
n2=facil_table_stat_high(2)*N2;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=facil_table_stat_low(3)*N1;
n2=facil_table_stat_high(3)*N2;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);


%[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3