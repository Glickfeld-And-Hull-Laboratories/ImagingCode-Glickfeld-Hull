

function MorphPlot(strEqnX, strEqnY)
%function MorphPlot(strEqnX, strEqnY)
%  Morphs all objects and axes limits
%
%strEqn should contain the variable 'x'
% Example  strEqn = 'log(2*x + 1)'


AxisF = {'Xlim' 'Ylim'   };
TickF = {'XTick' 'YTick' };
TickL = {'XTickLabel' 'YTickLabel'};
if nargin==1
    strEqn = {strEqnX  strEqnX};
elseif nargin==2
    strEqn = {strEqnX  strEqnY};
else
    error('Must supply an equation');
end

for f=1:length(AxisF)
    Axis{f} = get(gca, AxisF{f});
    Tick{f} = get(gca, TickF{f});
    Label{f} = get(gca, TickL{f});
end

DatF = {'XData' 'YData'};
H = get(gca,'Children');  %handles to children
for k=1:length(H)
    for f=1:length(DatF)
        try
            x = get(H(k), DatF{f});
            x = eval(strEqn{f});
            set(H(k), DatF{f}, x);
        catch
        end
    end
end

for f=1:length(AxisF)
    x = Axis{f};
    x = eval(strEqn{f});
    set(gca, AxisF{f}, x);
    
    x = Tick{f};
    x = eval(strEqn{f});
    set(gca, TickF{f}, x);
    
    %keep labels the same
    set(gca, TickL{f}, Label{f});
end
