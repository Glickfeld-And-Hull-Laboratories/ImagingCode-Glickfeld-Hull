%VB Gabor fitting, modified to fit 2D elliptical 

%params VB sends in: sta.kmean

s = struct;
s.data = data;
s.orig = b';
if P==2
s.dF = dF_mat(iCell,:);
s.ypos = i(iCell, :);
s.xpos = j(iCell, :);
end
m = @(pars,sftf) Gauss2D_ellipseMA(pars,sftf);
%s.m = func2str(m);

clear('x');
clear('y');
[m2,n2] = size(grid2.sfsf);
x(:,1) = log2(grid2.sfsf(:));
x(:,2) = log2(grid2.tftf(:));
Nperfreq = sqrt(size(x,1));
xhigh = reshape(x(:,1),Nperfreq,Nperfreq);
xhigh2 = interp2(xhigh,4);
Nperfreq2 = (size(xhigh2,1));
xhigh3 = reshape(xhigh2,Nperfreq2*Nperfreq2,1);

yhigh = reshape(x(:,2),Nperfreq,Nperfreq);
yhigh2 = interp2(yhigh,4);
Nperfreq2 = (size(yhigh2,1));
yhigh3 = reshape(yhigh2,Nperfreq2*Nperfreq2,1);
xhigh4 = [xhigh3 yhigh3];

uvar_high = zeros(Nperfreq2,2);
uvar_high(:,2) = xhigh2(:,1);
uvar_high(:,1) = yhigh2(1,:)';

x_plot = x;
y = s.data(:);

 Mask_NaN = isnan(y);
 ind_NaN = find(isnan(y));
 ind_noNaN = find(Mask_NaN == 0);
 if ~isempty(ind_NaN)
     x(ind_NaN,:) = [];
     y(ind_NaN) = [];
 end

sf_all = unique(x_plot(:,1));
tf_all = unique(x_plot(:,2));

%A sigma_SF sigma_TF sf0 tf0 xi
if P==1
    s.lb =  [.001 1 1  (min(x(:,1))) (min(x(:,2))) -2];
    s.ub =  [3 5   5    (max(x(:,1))) (max(x(:,2))) 2];
elseif P==2
s.lb = zeros(1,6);
s.ub = zeros(1,6);
s.lb(1,[1 2 3 6]) = [.001 .25 .25 -2];
s.ub(1,[1 2 3 6]) = [3 3 3 2];
if (max(x(:,1)))-(min(x(:,1)))>0
    s.lb(1,4) = min(x(:,1));
    s.ub(1,4) = max(x(:,1));
else
    sf_max = find(sf_all == max(x(:,1)));
    if sf_max == length(sf_all) 
        sf_max2 = sf_all(sf_max,:);
        sf_min2 = sf_all(sf_max-1,:);
    elseif sf_max == 1
        sf_max2 = sf_all(sf_max+1,:);
        sf_min2 = sf_all(sf_max,:);
    else
        sf_max2 = sf_all(sf_max+1,:);
        sf_min2 = sf_all(sf_max-1,:);
    end
    s.lb(1,4) = sf_min2;
    s.ub(1,4) = sf_max2;
end
if (max(x(:,2)))-(min(x(:,2)))>0
    s.lb(1,5) = min(x(:,2));
    s.ub(1,5) = max(x(:,2));
else
    tf_max = find(tf_all == max(x(:,2)));
    if tf_max == length(tf_all) 
        tf_max2 = tf_all(tf_max,:);
        tf_min2 = tf_all(tf_max-1,:);
    elseif tf_max == 1
        tf_max2 = tf_all(tf_max+1,:);
        tf_min2 = tf_all(tf_max,:);
    else
        tf_max2 = tf_all(tf_max+1,:);
        tf_min2 = tf_all(tf_max-1,:);
    end
    s.lb(1,5) = tf_min2;
    s.ub(1,5) = tf_max2;
end
end
%from before 110924 LG
% s.lb =  [.001 .25 .25  (min(x(:,1))) (min(x(:,2))) -2];
% s.ub =  [3 4   4    (max(x(:,1))) (max(x(:,2))) 2];

%changed on 111117 LG
% s.lb =  [.001 1 1  (min(x(:,1))) (min(x(:,2))) -2];
% s.ub =  [3 5   5    (max(x(:,1))) (max(x(:,2))) 2];
% 
%OLD WAS up to 110514:
% s.lb =  [.001 .25 .25  (min(x(:,1))-1) (min(x(:,2))-1) -2];
% s.ub =  [1000 3   3    (max(x(:,1))+1) (max(x(:,2))+1) 2];

%s.opt = 'silent';

%ph = 0:pi/2:3*pi/2;
%ori = 0:pi/4:3*pi/4;
%pos = [-15,-15;-15 0;-15 15;0 -15; 0 0;0 15; 15 -15; 15 0 ; 15 15];

if P==2
    Nsamps = 2;
elseif P==1
    Nsamps = 4; %changed from 2 to 4 on 111117 LG
end
dbin = [s.ub(1) - s.lb(1); ...
    s.ub(2) - s.lb(2); ...
    s.ub(3) - s.lb(3); ...
    s.ub(4) - s.lb(4); ...
    s.ub(5) - s.lb(5); ...
    s.ub(6) - s.lb(6)] ./ (Nsamps-1);
    
sigma_SF_vec = s.lb(2)+dbin(2)/2:dbin(2):s.ub(2);
sigma_TF_vec = s.lb(3)+dbin(3)/2:dbin(3):s.ub(3);
SF_vec = s.lb(4)+dbin(4)/2:dbin(4):s.ub(4);
TF_vec = s.lb(5)+dbin(5)/2:dbin(5):s.ub(5);
xi_vec = s.lb(6)+dbin(6)/2:dbin(6):s.ub(6);


clear temp;
index = 1;

% loop over initial parameters not to get stuck in local minimum

names = {'x2','resnorm'};
args = cell(5);

%options.Display = 'off';

for iSigSF = 1:length(sigma_SF_vec)
    for iSigTF = 1:length(sigma_TF_vec)
        for iSF = 1:length(SF_vec)
            for iTF = 1:length(TF_vec)
                for ixi = 1:length(xi_vec)
%                    fprintf('.');
                    s.x0 =  [max(max(s.data)) sigma_SF_vec(iSigSF) sigma_TF_vec(iSigTF) SF_vec(iSF) TF_vec(iTF) xi_vec(ixi)];

                    %                test = Gauss2D_ellipseMA(s.x0,x);

                    %                [args{:}] = lsqcurvefit(m,s.x0,grid2.sfsf,grid2.tftf,s.data,s.lb,s.ub,s.opt);
                    %                [args{:}] =
                    %                lsqcurvefit(m,s.x0,grid2.sfsf,grid2.tftf,s.data,s.lb,s.ub);
                    %                [args{:}] = lsqcurvefit(m,s.x0,x,y,s.lb,s.ub);
                    options = optimset('Display', 'off');
                    [x2,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB]  = lsqcurvefit(m,s.x0,x,y,s.lb,s.ub,options);
                    temp(index) = cell2struct({x2,Resnorm},names,2);
                    index = index +1;
                end
            end;
        end;
    end;
end


[val,ind]=min([temp.resnorm]);
tmp = [];
for countres = 1:length(temp)
    tmp = [tmp; temp(countres).resnorm temp(countres).x2];
end

s.fit = temp(ind);
s.x = s.fit.x2;
s.k2 = zeros(m2*n2,1);
s.k2(ind_noNaN) = m(s.fit.x2,x);
s.k2_plot = m(s.fit.x2,x_plot);

if P ==1
s.k2_plot_oversamp0 = m(s.fit.x2,xhigh4);
s.k2_plot_oversamp = reshape(s.k2_plot_oversamp0,Nperfreq2,Nperfreq2);
end

s.k2b = reshape(s.k2,m2,n2);
s.k2b_plot = reshape(s.k2_plot,m2,n2);

s.res = s.data-s.k2b;
s.Maxfit = max(max(s.k2b));
s.Maxdata = max(max(s.data));


%find highcut for SF and TF 
%x(:,1) = log2(grid2.sfsf(:));
%x(:,2) = log2(grid2.tftf(:));
%    TF_vec0 = ind_TFuse(:,2);
%    SF_vec0 = ind_SFuse(:,2);
x00 = zeros(size(grid2.sfsf00,1)*size(grid2.sfsf00,2),2);
x00(:,1) = (grid2.sfsf00(:));
x00(:,2) = (grid2.tftf00(:));
k2b00 = m(s.fit.x2,x00);
%now find 10% and 50%
s.Maxfit00 = max(max(k2b00));

%tmp2 = max(tmp,[],1); %max over TFs
%tmp3 = tmp > .5*max(max(tmp));
%imagesc(2.^TF_vec00,flipud(2.^SF_vec00),flipud(tmp3'))
MaxSF00 = s.x(4);
MaxTF00 = s.x(5);

indSF50 = find(k2b00>.5*s.Maxfit00 & x00(:,1)>MaxSF00);
indTF50 = find(k2b00>.5*s.Maxfit00 & x00(:,2)>MaxTF00);
SFhicut_50 = NaN;
TFhicut_50 = NaN;
if ~isempty(indSF50)
    SFhicut_50 = max(x00(indSF50,1));
end
if ~isempty(indTF50)
    TFhicut_50 = max(x00(indTF50,2));
end
s.SFhicut_50 = SFhicut_50;
s.TFhicut_50 = TFhicut_50;

%repeat for 10% cutoff

indSF10 = find(k2b00>.1*s.Maxfit00 & x00(:,1)>MaxSF00);
indTF10 = find(k2b00>.1*s.Maxfit00 & x00(:,2)>MaxTF00);
SFhicut_10 = NaN;
TFhicut_10 = NaN;
if ~isempty(indSF10)
    SFhicut_10 = max(x00(indSF10,1));
end
if ~isempty(indTF10)
    TFhicut_10 = max(x00(indTF10,2));
end
s.SFhicut_10 = SFhicut_10;
s.TFhicut_10 = TFhicut_10;




if PLOTIT_FIT == 1
if  P ==1
    h = subplot(sum(H_ttest(:,iRun),1),4,start);
    text(0.5, 1, area_list(iCell,:));
    xlim([0 2])
    ylim([0 2])
    axis off
    h = subplot(sum(H_ttest(:,iRun),1),4,1+start);
    imagesc(s.data); colormap('gray');caxis([0 max(max(Im_mat_USE, [],1),[],2)]); axis image;
    set(h,'XTick',[1:size(uvar,1)])
    set(h,'YTick',[1:size(uvar,1)])
	set(h,'YTickLabel',flipud(uvar(:,2)));set(h,'XTickLabel',(uvar(:,1)))
    h = subplot(sum(H_ttest(:,iRun),1),4, 2+start);
    imagesc(s.k2b_plot); colormap('gray'); caxis([0 max(max(Im_mat_USE, [],1),[],2)]); axis image;
    set(h,'XTick',[1:size(uvar,1)])
    set(h,'YTick',[1:size(uvar,1)])
    set(h,'YTickLabel',flipud(uvar(:,2)));set(h,'XTickLabel',(uvar(:,1)))
    h = subplot(sum(H_ttest(:,iRun),1),4, 3+start);
    imagesc(s.k2_plot_oversamp); colormap('gray'); caxis([0 max(max(Im_mat_USE, [],1),[],2)]); axis image;
    set(h,'XTick',[1:size(uvar_high,1)])
    set(h,'YTick',[1:size(uvar_high,1)])
    set(h,'YTickLabel',flipud(uvar_high(:,2)));set(h,'XTickLabel',(uvar_high(:,1)));
%     xtxt = interp1(str2num(get(h,'XTickLabel')),get(h,'XTick')', s.x(5));
%     ytxt = interp1(str2num(get(h,'YTickLabel')),flipud(get(h,'YTick')'), s.x(4));
%     hold on
% %     plot(xtxt,ytxt,'k.');
    axis off
    start = start+4;
if iCell == nCells
    suptitle([mouse ' ' date ' ' str_run(iRun,:)])
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_' str_run(iRun,:) '_SFxTF_fits.ps']);
	print(gcf, '-depsc', fn_out);
end
end

if P==2
    if start >100
        figure;
        start = 1;
    end
    h = subplot(10,10, start);
    MIN = 0;
    MAX = max(max(Im_mat_USE(iCell,:),[],2),[],1);
    imagesq(s.data,[MIN MAX]); colormap('gray'); axis image;
    colormap(gray)
    h = subplot(10,10, start+1);
    MIN = 0;
    MAX = max(max(Im_mat_USE(iCell,:),[],2),[],1);
    imagesq(s.k2b_plot,[MIN MAX]); colormap('gray'); axis image;
    colormap(gray)
    start = start+2;
end

end

if SAVEALLDATA == 0
    s = rmfield(s,'res');
    s = rmfield(s,'k2b');
    s = rmfield(s,'k2b_plot');
    s = rmfield(s,'k2_plot');
    s = rmfield(s,'k2');
    s = rmfield(s,'data');
    s = rmfield(s,'Maxfit00');
    s = rmfield(s,'Maxfit');
    s = rmfield(s,'Maxdata');
    s = rmfield(s,'orig');
    s = rmfield(s,'ub');
    s = rmfield(s,'lb');
    s = rmfield(s,'fit');
    s = rmfield(s,'dF');
    s = rmfield(s,'xpos');
    s = rmfield(s,'ypos');
end



%this.gab1 = s;


