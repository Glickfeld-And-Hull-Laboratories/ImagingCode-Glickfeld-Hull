
function [xm,xs, ym,ys] = PlotBinnedData2(Bins,X,Y,N, SortWeightXY, clr, sym, lw, bUseMedian, bPlotRaw)
%function [xm,xs,xs,ym,ys,ys] = PlotBinnedData(Bins,X,Y,N, SortWeightXY, clr, sym, lw, bUseMedian, bPlotRaw)
%
%Bins = two options
%   scalar = number of data bins
%   vect   = bin boundaries based on SortWeighted data (note, data smaller than first and larger than last boundary will be discarded; use [-inf x x x +inf] to ensure all data is kept)
%
%X,Y = vector X/Y data
%
%N (optional) = fractional statistical "n" for each data pair.  Set to [] to default to n=1 for each data pair.
%
%SortWeightXY.  (optional)
%  To sort only on X = [1 0]  (default)
%  To sort only on Y = [0 1]
%  To sort on both recommend = [median(X)  median(Y)].
%  Note: can set either to negative to invert order
%
%bUseMedian = 1 to use median (default) or 0 to use mean

%default params if necessary
try
    bPlotRaw;
catch
    bPlotRaw = 1;
end
try
    bUseMedian;
catch
    bUseMedian = 1;
end
try
    N;
catch
    N = [];
end
try
    SortWeightXY;
    if isempty(SortWeightXY)
        error;
    end
catch
    SortWeightXY = [1 0];
end
try
    sym;
catch
    sym = 'o';
end
try
    clr;
catch
    clr = [0 0 0];
end
try
    lw;
catch
    lw = 2;
end


%error-check input data
X = X(:);
Y = Y(:);
N = N(:);
if isempty(N)
    N = ones(length(X),1);
end
if length(X)~=length(Y) || length(X)~=length(N)
    error('X and Y (and N if not empty) must be vectors of the same length');
end


%sort
SortWeightXY = SortWeightXY/sum(SortWeightXY); %normalize
[SW,I] = sort( X*SortWeightXY(1) + Y*SortWeightXY(2) );
X = X(I);
Y = Y(I);
N = N(I);

if length(Bins)==1
    for q=1:Bins
        J{q} = find(cumsum(N)>sum(N)*(q-1)/Bins & cumsum(N)<=sum(N)*q/Bins);
        if bUseMedian
            [~, xs(q), ~, xm(q)] = WeightedMeanSEM(X(J{q}),N(J{q}));
            [~, ys(q), ~, ym(q)] = WeightedMeanSEM(Y(J{q}),N(J{q}));
        else
            [xm(q), xs(q)] = WeightedMeanSEM(X(J{q}),N(J{q}));
            [ym(q), ys(q)] = WeightedMeanSEM(Y(J{q}),N(J{q}));
        end
    end
else
    for q=1:length(Bins)-1
        J{q} = find(SW>Bins(q) & SW<=Bins(q+1));
        if bUseMedian
            [~, xs(q), ~, xm(q)] = WeightedMeanSEM(X(J{q}),N(J{q}));
            [~, ys(q), ~, ym(q)] = WeightedMeanSEM(Y(J{q}),N(J{q}));
        else
            [xm(q), xs(q)] = WeightedMeanSEM(X(J{q}),N(J{q}));
            [ym(q), ys(q)] = WeightedMeanSEM(Y(J{q}),N(J{q}));
        end
    end
end
xscale = get(gca,'xscale'); yscale = get(gca,'yscale');  %see what it is
set(gca,'xscale',xscale,'yscale',yscale);  %ensure it stays the same
ax = axis;
if strcmpi(get(gca,'xscale'),'log')
    dX = -1.2; %(ax(2)/ax(1)).^(1/10);  %log
else
    dX = (ax(2)-ax(1))/75;  %linear
end
if strcmpi(get(gca,'yscale'),'log')
    dY = -1.2;  %log
else
    dY = (ax(4)-ax(3))/75;
end
if bPlotRaw==2
    %color-mapping
    for w=1:2
        for q=1:length(J)
            if length(J)==4
                clr = MapValue2Color(q,[1.5 length(J)-0.5], jet);
            elseif length(J)==3
                clr = MapValue2Color(q,[1.6 length(J)], jet);
            else
                clr = MapValue2Color(q,[1 length(J)], jet);
            end
            %clr = [0 0 0]+(mod(q,2))*0.5;
            wclr  = clr*0.5 + [1 1 1]*0.5;
            wclr2 = clr*0.2 + [1 1 1]*0.8;
            if w==1
                %scatter first
                plot(X(J{q}),Y(J{q}),sym,'color',wclr, 'linestyle', 'none', 'markerfacecolor', wclr2); hold on;
            else
                %binned data second
                if q>1
                    plot(xm(q-1:q),ym(q-1:q),'-','color',clr,'MarkerSize',15,'linewidth',lw); hold on;
                end
                plot(xm(q),ym(q),sym,'color',clr,'MarkerSize',15,'linewidth',lw); hold on;
                ErrorPlot(gca,xm(q),xs(q),xs(q),ym(q),ys(q),ys(q), dX,dY,lw, clr);
            end
        end
    end
else
    %normal mode
    if bPlotRaw==1
        wclr = clr;
        try
            if ~ischar(clr)
                wclr = wclr*0.4 + [1 1 1]*0.6;
            end
        catch
        end
        plot(X,Y,sym,'color',wclr, 'linestyle', 'none'); hold on;
    end
    plot(xm,ym,sym,'color',clr,'MarkerSize',15,'linewidth',lw); hold on;
    ErrorPlot(gca,xm,xs,xs,ym,ys,ys, dX,dY,lw,clr);
end


