% Size tuning fitting by KM, modified from Fit_2Dellipse_LG_Ret_KM

% start with dumdum (dF/F responses for single cell+contrast)
% and szs0 (stim size for each response)

s = struct;
s.data = dumdum;
s.szs0 = szs0;

m1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))));
m2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))));
s.m1 = func2str(m1);
s.m2 = func2str(m2);

% lower/upper bounds
%Ae ke xe
s.lb1 =  [0 0 0];
s.ub1 =  [inf 1 80];
%Ae k1 x1 Ai k2 x2
s.lb2 =  [0 0 0 0 0 0];
s.ub2 =  [inf 1 80 inf 1 80];

% for sampling different initial guesses
Nsamps = 2; %changed from 2 to 3
dbin1 = [s.ub1(1) - s.lb1(1); ...
    s.ub1(2) - s.lb1(2); ...
    s.ub1(3) - s.lb1(3)]./ (Nsamps-1);
dbin2 = [s.ub2(1) - s.lb2(1); ...
    s.ub2(2) - s.lb2(2); ...
    s.ub2(3) - s.lb2(3); ...
    s.ub2(4) - s.lb2(4); ...
    s.ub2(5) - s.lb2(5); ...
    s.ub2(6) - s.lb2(6)] ./ (Nsamps-1);

ke_vec = s.lb1(2)+dbin1(2)/2:dbin1(2):s.ub1(2);
xe_vec = s.lb1(3)+dbin1(3)/2:dbin1(3):s.ub1(3);

k1_vec = s.lb2(2)+dbin2(2)/2:dbin2(2):s.ub2(2);
x1_vec = s.lb2(3)+dbin2(3)/2:dbin2(3):s.ub2(3);
k2_vec = s.lb2(5)+dbin2(5)/2:dbin2(5):s.ub2(5);
x2_vec = s.lb2(6)+dbin2(6)/2:dbin2(6):s.ub2(6);

% loop over initial parameters not to get stuck in local minimum
% start with model1
clear temp;
index = 1;
names = {'c1','OF1','x0'};
% model1
for ike = 1:length(ke_vec)
    for ixe = 1:length(xe_vec)
        % fprintf('.');
        guess1 = [0.5*maxMean 0.35 x_max];
        %guess1 = [0.5*maxMean ke_vec(ike) xe_vec(ixe)];
        %guess1 = [0.5*maxMean 0.35 xe_vec(ixe)];
        
        [c1,OF1] = fminsearchbnd(@(c)sigevalpen1(c,s.szs0,s.data),guess1,s.lb1,s.ub1,opts);
        temp(index) = cell2struct({c1,OF1,guess1},names,2);
        index = index +1;
    end
end
[val,ind]=min([temp.OF1]);
s.fit1 = temp(ind);

% model2
clear temp;
index = 1;
names = {'c2','OF2','x0'};
for ik1 = 1:length(k1_vec)
    for ix1 = 1:length(x1_vec)
        for ik2 = 1:length(k2_vec)
            for ix2 = 1:length(x2_vec)
                % fprintf('.');
                %guess2 = [maxMean k1_vec(ik1) x1_vec(ix1) maxMean-mean(s.data) k2_vec(ik2) x2_vec(ix2)];
                guess2 = [maxMean 0.4 x_max/2 maxMean-mean(dumdum) 0.3 1.5*x_max];
                %guess2 = [maxMean 0.4 x1_vec(ix1) maxMean-mean(dumdum) 0.3 x2_vec(ix2)];
                
                [c2,OF2] = fminsearchbnd(@(c)sigevalpen2(c,s.szs0,s.data),guess2,s.lb2,s.ub2,opts);
                temp(index) = cell2struct({c2,OF2,guess2},names,2);
                index = index +1;
            end
        end
    end
end
[val,ind]=min([temp.OF2]);
s.fit2 = temp(ind);

% with selected fits, compute statistics
s.res1 = s.data-m1(s.fit1.c1,s.szs0);
s.res2 = s.data-m2(s.fit2.c2,s.szs0);
s.SStot = sum((s.data-mean(s.data)).^2);
s.SSE1 = norm(s.res1).^2; % SSE is squared 2-norm of error
s.SSE2 = norm(s.res2).^2;
s.Rsq1 = 1 - s.SSE1/s.SStot;
s.Rsq2 = 1 - s.SSE2/s.SStot;

szRng = linspace(0,max(szs)); % oversample at default 100 points
% find prefSize/SI for model1 (pref at 90% max)
s.fitout1 = m1(s.fit1.c1,szRng);
s.maxResp1 = max(s.fitout1);
s.prefSize1 = szRng(find(s.fitout1>(0.9*s.maxResp1),1));
if ~length(s.prefSize1)
    s.prefSize1=max(szs);
end
s.suppInd1 = 0;

% find prefSize/SI for model2 (pref at max)
s.fitout2 = m2(s.fit2.c2,szRng);
s.maxResp2 = max(s.fitout2);
s.prefSize2 = szRng(find(s.fitout2==s.maxResp2,1));
s.suppInd2 = 1 - s.fitout2(end)/s.maxResp2;

% F-test
dfS = length(s.data) - 3; % single sigmoid model df = #trials - 3 model parameters
dfD = length(s.data) - 6; % double sigmoid model df = #trials - 6 model parameters
s.Fcrit = finv(0.95,dfD,dfS); % measure critical F at a=0.05
s.Fscore = (s.SSE2 - s.SSE1)/(dfD-dfS) ./ (s.SSE2/dfD);
s.Ftest = s.Fscore>s.Fcrit;
% check if Ai or ki below threshold, if so revert to m1
% Ai < 1% of maxResp, or ki <1e-4
if (s.fit2.c2(4)<0.01*s.maxResp2) || (s.fit2.c2(5)<1e-4)
    s.Ftest = 0;
end

% store values from selected test (0=single,1=double)
if s.Ftest
    s.prefSize = s.prefSize2;
    s.suppInd = s.suppInd2;
    s.Fstr1=''; s.Fstr2='*';
else
    s.prefSize = s.prefSize1;
    s.suppInd = s.suppInd1;
    s.Fstr1='*'; s.Fstr2='';
end

% plot both fits (with initial guesses) and residual (4 plots each, in blocks)
if PLOTIT_FIT == 1
    if start ==37
        set(gcf, 'Position', [0 0 800 1000]);
        fn_out = fullfile(fnout, datemouse, datemouserefrun, [datemouserefrun '_SizeTuneFits' num2str(ifig) '.pdf']);
        print(fn_out,'-dpdf')
        figure;
        ifig = 1+ifig;
        start = 1;
    end
    h = subplot(6,6,start);
    if nCon>1
        errorbar([0 szs],[0 sizeMean(:,nCon,iCell)'],[0 sizeSEM(:,nCon,iCell)'])
    elseif nrun>1
        errorbar([0 szs],[0 sizeMean(:,1,iCell)'],[0 sizeSEM(:,1,iCell)'])
    end
    hold on
    plot(s.szs0,s.data,'.b')
    if s.Ftest
        plot(szRng,s.fitout1,':k')
        plot(szRng,s.fitout2,'-r')
    else
        plot(szRng,s.fitout1,'-r')
        plot(szRng,s.fitout2,':k')
    end
    plot(szRng,m1(s.fit1.x0,szRng),'g--')
    plot(szRng,m2(s.fit2.x0,szRng),'g--')
    hold off
    ylim([min([-0.5*s.maxResp1 min(s.data)]) 1.2*max([s.maxResp2 max(s.data)])])
    title([num2str(iCell), ' dist:' num2str(cellDists(iCell),3) s.Fstr2]);
    
    % residual plots
%     plot(s.szs0,s.res1,'.')
%     title(['R^2_1 ' num2str(s.Rsq1,3)])
%     plot(s.szs0,s.res2,'.')
%     title(['R^2_2 ' num2str(s.Rsq2,3)])
    
    start = start+1;
end

if SAVEALLDATA == 0
    s = rmfield(s,'data');
    s = rmfield(s,'szs0');
    s = rmfield(s,'lb1');
    s = rmfield(s,'ub1');
    s = rmfield(s,'lb2');
    s = rmfield(s,'ub2');
    s = rmfield(s,'fitout1');
    s = rmfield(s,'fitout2');
    s = rmfield(s,'res1');
    s = rmfield(s,'res2');
end
