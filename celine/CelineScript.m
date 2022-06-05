function CelineScript(Qvect)
%function CelineScript(QTransformValues)
%
%Scatter plot, with data binned.
%Note: nonlinear transform is possible --> log10(Q*DATA    + 1);

%  Q ranges from 0 to 8 or so.   The larger, the more log-like it becomes.
%
%Input: 
%  Qvect;  set of Q values to plot (animated)
%larger = more spread

load DART_dFoF_data

if nargin == 0
    fprintf('Running Video mode;\n  Try typing: >> CelineScript(0); for linear.\n  Try typing: >> CelineScript(9); for most log-like\n\n');
    Qvect = [0:0.5:9   9:-0.5:0];
end

set(0,'DefaultFigureWindowStyle','normal');

for Q = Qvect  
    bShowColors = 1; %set to 1 for colors
    ZoomOnHTP = 1;
        
    DAT=stat_resp;
    
    DATmin = min(DAT);
    DATmax = max(DAT);
    LIM = [...
        min(DATmin(2), DATmin(1)/ZoomOnHTP) ...
        max(DATmax(2), DATmax(1)/ZoomOnHTP)];
    LIM = { LIM,     LIM/ZoomOnHTP};   %axis limits
    TICK = -0.1:0.1:0.6;   TICK= { TICK          TICK/ZoomOnHTP};
    Title = {'-HTP'                           '+HTP'};
    
    %bin boundaries
    BIN=    {[-inf  0.01  0.1  0.19  inf],    [-inf  0.01   0.05  0.11 inf]};
    
    
    clf;
    for k=1:2  %HT- then HT+
        subplot(1,2,k); cla; hold on;
        %formating code
        
        DATA = stat_resp;
        set(gca,'xtick',TICK{k});
        set(gca,'ytick',TICK{k});
        title([Title{k}]);
        axis equal; axis([LIM{k} LIM{k}]);
        xlabel('baseline activity');
        ylabel('YM90K.1^{DART.2}');
        
        set(gca,'fontsize',12);
        grid on;
        
        X = DATA(find(HT_ind==k-1),2);   %before YM90K
        Y = DATA(find(HT_ind==k-1),1);    %after YM90K
        
        %raw data
        %plot(X,Y,'o', 'color', [0.6 0.6 1]);  hold on;
        
        %reference lines
        tmp=linspace(LIM{k}(1), LIM{k}(2), 20);
        plot(tmp,tmp,'k-', 'linewidth', 1); hold on; %y=x line
        
        %binned data
        PlotBinnedData2(BIN{k}, X, Y,[],[std(X)  std(Y)],[0 0 1],'-o',2,0,1+bShowColors);
        PlotBinnedData2(1, X, Y,[],[],'k','-',2,0,0);
        
        if Q>0
            MorphPlot(['log(' num2str(Q) '*x + 1)']);
            title([Title{k} ',  = log(' num2str(Q) 'x + 1)']);
            ax=axis; ax([1 3])=max(ax(1),-0.5); axis(ax);
        end
    end
    drawnow;
end
