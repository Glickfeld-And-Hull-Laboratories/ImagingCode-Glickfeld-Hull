% Fit_2Dellipse_ret_CC.m
% 2D elliptical Gaussian fit for retinotopic mapping
% Adapted from Fit_2Dellipse_LG_Ret.m
%
% Parameters (s.x):
%   (1) A         - amplitude
%   (2) sigma_El  - elevation width (deg)
%   (3) sigma_Az  - azimuth width (deg)
%   (4) El0       - elevation center (deg)
%   (5) Az0       - azimuth center (deg)
%   (6) xi        - rotation (fixed at 0)
%
% Expects in workspace: data, b, grid2, isilonName, mouse, date,
%                       ret_run_str, PLOTIT_FIT, SAVEALLDATA,
%                       h_all, iCell, ifig, start

s = struct;
s.data = data;
s.orig = reshape(b', size(b,1)*size(b,2), 1);

m = @(pars,sftf) Gauss2D_ellipseMA(pars,sftf);

clear x y
[m2,n2] = size(grid2.AzAz);

x = [grid2.AzAz(:) grid2.ElEl(:)];
uvar.Az = unique(x(:,1));
uvar.El = unique(x(:,2));
xNperfreq = size(grid2.AzAz00,1);
yNperfreq = size(grid2.AzAz00,2);
xyhigh = [reshape(grid2.AzAz00, xNperfreq*yNperfreq, 1), ...
          reshape(grid2.ElEl00, xNperfreq*yNperfreq, 1)];

x_plot = x;
y = s.data(:);

CM(1) = sum(x(:,1).*s.orig) / sum(s.orig);
CM(2) = sum(x(:,2).*s.orig) / sum(s.orig);

ind_NaN   = find(isnan(y));
ind_noNaN = find(~isnan(y));
if ~isempty(ind_NaN)
    x(ind_NaN,:) = [];
    y(ind_NaN)   = [];
end

% bounds: A  sig_El  sig_Az  El0           Az0           xi
s.lb = [.001   1       1     min(x(:,1))   min(x(:,2))   0];
s.ub = [3    100     100     max(x(:,1))   max(x(:,2))   0];

Nsamps = 2;
dbin = (s.ub - s.lb)' ./ (Nsamps-1);

sigma_El_vec = s.lb(2) + dbin(2)/2 : dbin(2) : s.ub(2);
sigma_Az_vec = s.lb(3) + dbin(3)/2 : dbin(3) : s.ub(3);
El_vec       = s.lb(4) + dbin(4)/2 : dbin(4) : s.ub(4);
Az_vec       = s.lb(5) + dbin(5)/2 : dbin(5) : s.ub(5);

clear temp
index = 1;
names = {'x2','resnorm'};
options = optimset('Display','off');

for iSigEl = 1:length(sigma_El_vec)
    for iSigAz = 1:length(sigma_Az_vec)
        for iEl = 1:length(El_vec)
            for iAz = 1:length(Az_vec)
                s.x0 = [max(s.data(:)) sigma_El_vec(iSigEl) sigma_Az_vec(iSigAz) El_vec(iEl) Az_vec(iAz) 0];
                [x2, Resnorm] = lsqcurvefit(m, s.x0, x, y, s.lb, s.ub, options);
                temp(index) = cell2struct({x2, Resnorm}, names, 2);
                index = index + 1;
            end
        end
    end
end

[~, ind] = min([temp.resnorm]);
s.fit = temp(ind);
s.x   = s.fit.x2;

s.k2              = zeros(m2*n2, 1);
s.k2(ind_noNaN)   = m(s.x, x);
s.k2_plot         = m(s.x, x_plot);
s.k2_plot_oversamp0 = m(s.x, xyhigh);
s.k2_plot_oversamp  = reshape(s.k2_plot_oversamp0, xNperfreq, yNperfreq);
s.k2b             = reshape(s.k2, m2, n2);
s.k2b_plot        = reshape(s.k2_plot, m2, n2);
s.res             = s.data - s.k2b;
s.Maxfit          = max(s.k2b(:));
s.Maxdata         = max(s.data(:));

% high-cut on oversampled grid
x00(:,1) = grid2.ElEl00(:);
x00(:,2) = grid2.AzAz00(:);
k2b00        = m(s.x, x00);
s.Maxfit00   = max(k2b00);
MaxEl00      = s.x(4);
MaxAz00      = s.x(5);

indEl50 = find(k2b00 > .5*s.Maxfit00 & x00(:,1) > MaxEl00);
indAz50 = find(k2b00 > .5*s.Maxfit00 & x00(:,2) > MaxAz00);
s.Elhicut_50 = NaN; s.Azhicut_50 = NaN;
if ~isempty(indEl50); s.Elhicut_50 = max(x00(indEl50,1)); end
if ~isempty(indAz50); s.Azhicut_50 = max(x00(indAz50,2)); end

indEl10 = find(k2b00 > .1*s.Maxfit00 & x00(:,1) > MaxEl00);
indAz10 = find(k2b00 > .1*s.Maxfit00 & x00(:,2) > MaxAz00);
s.Elhicut_10 = NaN; s.Azhicut_10 = NaN;
if ~isempty(indEl10); s.Elhicut_10 = max(x00(indEl10,1)); end
if ~isempty(indAz10); s.Azhicut_10 = max(x00(indAz10,2)); end

if PLOTIT_FIT == 1
    if start == 65
        fn_out = fullfile(isilonName, '/home/ACh/Analysis/2P_analysis', mouse, date, RetImgFolder, ...
            [date '_' mouse '_' ret_run_str '_RFfits' num2str(ifig) '.pdf']);
        print(fn_out, '-dpdf')
        figure;
        ifig  = ifig + 1;
        start = 1;
    end
    MAX  = max(data(:));
    dF   = num2str(MAX*100);
    if h_all(1,iCell); sig_str = ' **'; else; sig_str = []; end

    h = subplot(8,8, start);
    imagesc(s.data); colormap('gray'); axis image;
    colorbar('YTick', MAX, 'YTickLabel', dF(1:3))
    title([num2str(iCell) sig_str]);
    set(h, 'XTick', 1:size(uvar.Az,1), 'YTick', 1:size(uvar.El,1), ...
        'YTickLabel', flipud(uvar.El), 'XTickLabel', uvar.Az)

    h = subplot(8,8, start+1);
    imagesc(s.k2b_plot); colormap('gray'); axis image; axis off;
    set(h, 'XTick', 1:size(uvar.Az,1), 'YTick', 1:size(uvar.El,1), ...
        'YTickLabel', flipud(uvar.El), 'XTickLabel', uvar.Az)

    h = subplot(8,8, start+2);
    imagesc(s.k2_plot_oversamp); colormap('gray'); axis image; axis off;

    h = subplot(8,8, start+3);
    imagesc(s.res); colormap('gray'); axis image; axis off;
    set(h, 'XTick', 1:size(uvar.Az,1), 'YTick', 1:size(uvar.El,1), ...
        'YTickLabel', flipud(uvar.El), 'XTickLabel', uvar.Az)
    clim([0 MAX])

    start = start + 4;
end

if SAVEALLDATA == 0
    s = rmfield(s, {'res','k2b','k2','data'});
end