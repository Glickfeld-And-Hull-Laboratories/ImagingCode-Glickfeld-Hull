

function ErrorPlot(ax, X, XL, XU, Y, YL, YU, BarWidth, BarHeight, Linewidth, Color)
%function ErrorPlot(ax, X, XL, XU, Y, YL, YU, BarWidth, BarHeight, Linewidth, Color)
%
%NOTE: 
%if data is plotted on a linear scale,
%    then set BarWidth/BarHeight to positive values (additive +-W; +-H)
%if data is plotted on a log scale, 
%    then set BarWidth/BarHeight to negative values (mult  */W; */H)

X = X(:);
if isempty(XU)
    XU = X*0;
elseif length(XU)==1
    XU = X*0 + XU;
else
    XU = XU(:);
end

if isempty(XL)
    XL = X*0;
elseif length(XL)==1
    XL = X*0 + XL;
else
    XL = XL(:);
end
Y = Y(:);
if isempty(YU)
    YU = Y*0;
elseif length(YU)==1
    YU = Y*0 + YU;
else
    YU = YU(:);
end
if isempty(YL)
    YL = Y*0;
elseif length(YL)==1
    YL = Y*0 + YL;
else
    YL = YL(:);
end
BarWidth = BarWidth(:);
BarHeight = BarHeight(:);

tmp =  [length(X) length(XL) length(XU) length(Y) length(YL) length(YU)];
if min(tmp) ~= max(tmp)
    error('all vectors must be the same length');
end

N = length(X); %number of points
XL = X-XL;  %set XL/XU to absolute values (not just deviation from X)
XU = X+XU;
YL = Y-YL;
YU = Y+YU;
if length(BarHeight)==1
    BarHeight = Y*0 + BarHeight;
end
if any(BarHeight>=0)
    %linear addition
    HU = Y*0 + BarHeight/2;
    HL = HU;
else
    %multiplicative factors
    BarHeight = -BarHeight; %remove negative flag
    HU = Y.*sqrt(BarHeight) - Y;
    HL = Y - Y./sqrt(BarHeight);
end

if length(BarWidth)==1
    BarWidth = X*0 + BarWidth;
end
if any(BarWidth>=0)
    %linear addition
    WU = X*0 + BarWidth/2;
    WL = WU;
else
    %multiplicative factors
    BarWidth = -BarWidth; %remove negative flag
    WU = X.*sqrt(BarWidth) - X;
    WL = X - X./sqrt(BarWidth);
end


%each data point gets six line segments, 
%each line segment is represented by 3 x/y points 
%  (ends of line segement, separated by NaN's)
xx = zeros(18*N,1);
yy = zeros(18*N,1);
for k=1:N
   xx(18*(k-1)+1 : 18*k) = [XL(k) XU(k) nan   X(k)  X(k) nan  XL(k)  XL(k)  nan  XU(k)  XU(k)  nan  X(k)-WL(k) X(k)+WU(k) nan  X(k)-WL(k) X(k)+WU(k) nan];
   yy(18*(k-1)+1 : 18*k) = [ Y(k)  Y(k) nan  YL(k) YU(k) nan  Y(k)-HL(k) Y(k)+HU(k) nan  Y(k)-HL(k) Y(k)+HU(k) nan  YL(k)  YL(k)  nan  YU(k)  YU(k)  nan];
end

try
    plot(ax, xx, yy, 'linewidth', Linewidth, 'color', Color);
catch
    'debug'
end
