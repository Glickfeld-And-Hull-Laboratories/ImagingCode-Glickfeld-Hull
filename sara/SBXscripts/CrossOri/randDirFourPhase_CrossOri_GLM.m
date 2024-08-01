close all;clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary')
doPlot = 1;
ds = ['CrossOriRandDirFourPhase_ExptList_SG'];
svName = 'randPhase';
eval(ds)
driver = 'SLC';
img_area = {'V1';'L2/3'}; %LM
inj_area = 'V1';

%%
iexp = [6];
mouse = expt(iexp).mouse;
date = expt(iexp).date;
if isfield(expt,'copFolder') 
    ImgFolder = expt(iexp).copFolder;
else
    ImgFolder = expt(iexp).coFolder;
end
nrun = length(ImgFolder);
run_str = 'runs-002';

fprintf([mouse ' ' date '\n'])

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_varsforGLM.mat']))
nCells = size(resp_cell{1,1,1},1);
nDir = size(resp_cell,1);
ntrials = size(stimDir_all,2);
dirs = unique(stimDir_all);

%%
respGrat = resp_cell(:,1,1);
trialsGrat = trialInd(:,1,1);

resp_list = [];
trial_list = [];
for id = 1:nDir
    resp_list = [resp_list, respGrat{id}];
    trial_list = [trial_list, trialsGrat{id}];
end
dir_list = stimDir_all(trial_list);

for id = 1:nDir
    dir_list(find(dir_list==dirs(id)))=id;
end
resp_list = resp_list';
dir_list = dir_list';





save(fullfile(base, 'Analysis\2P\CrossOri\test', [date '_' mouse '_' run_str '_glmData.mat']), 'date','mouse','dir_list','resp_list');




%% Multinominal logistic regression
% B = mnrfit(X,Y)
% [B,dev,stats] = mnrfit(___)
% 
% B = mnrfit(X,Y) returns a matrix, B, of coefficient estimates for a multinomial logistic regression of 
% the nominal responses in Y on the predictors in X.

% Input Arguments:
%  -- X: Observations on predictor variables, specified as an n-by-p matrix. X contains n observations for 
%     p predictors.
%  -- Y: Response category labels, specified as an n-by-k matrix, or an n-by-1 numeric vector, string vector, 
%     categorical array, or cell array of character vectors.
resp_list = resp_list(:,1:4);

[B, dev, stats] = mnrfit(resp_list,x,'Interactions','off');

%% Support vector machine
% x = randperm(450,450)';
% train_resp_list = resp_list(x(200),:);
% test_resp_list = resp_list(x(201:450),:);
% train_dir_list = dir_list(x(1:200));
% test_dir_list = dir_list(x(201:450));

Mdl = fitcecoc(resp_list,dir_list,'Holdout',0.1);

[label, score] = predict(Mdl,test_resp_list);

idx = [];
for i = 1:length(label)
    idx(i) = label(i) == test_dir_list(i);
end
perCorr = length(find(idx==1))/length(idx)


%% GLM fit
x = randperm(450,450)';
train_resp_list = resp_list(x(1:400),:);
test_resp_list = resp_list(x(401:450),:);
train_dir_list = dir_list(x(1:400));
test_dir_list = dir_list(x(401:450));

x = train_resp_list;
y = train_dir_list;
xpred = test_resp_list;

% y(find(y==12)) = 0;
[b, dev, stats ]= mnrfit(x, y);


[b, dev, stats ]= glmfit(x, y,'binomial','link','logit');

yhat = glmval(b,xpred,'logit');

%% GLM fit
c = find(dir_list == 3);
d = find(dir_list == 7);
e = find(dir_list == 12);
x = [c(1:20);d(1:20);e(1:20)];
t = [c(21:35);d(21:35);e(21:35)];
train_resp_list = resp_list(x,1:10);
test_resp_list = resp_list(t,1:10);
train_dir_list = dir_list(x);
test_dir_list = dir_list(t);


x = train_resp_list;
y = train_dir_list;
y(find(y==3))=1;
y(find(y==7))=2;
y(find(y==12))=3;
xpred = test_resp_list;

dY = dummyvar(y);



[beta1, dev1, stats1] = glmfit(x, dY, 'binomial', 'link', 'logit');


[B, dev, stats] = mnrfit(x,dY,'model','ordinal','Interactions','off');





%testing
yfit = glmval(b,test_resp_list,'probit');

