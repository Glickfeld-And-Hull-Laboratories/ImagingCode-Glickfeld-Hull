function setFigParams4Print(paperOriString)
%% set figures to resize to 8.5x11 for printing to pdf, in either landscape or portrait orientation, as indicated by paperOriString
set(gca,'fontsize',8)
if strcmp(paperOriString,'portrait')
    set(0,'defaultfigurepaperorientation','portrait');
    set(0,'defaultfigurepapersize',[8.5 11]);
    set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
elseif strcmp(paperOriString,'landscape')
    set(0,'defaultfigurepaperorientation','landscape');
    set(0,'defaultfigurepapersize',[11 8.5]);
    set(0,'defaultfigurepaperposition',[.25 .25 ([11 8.5])-0.5]);    
else
    error('input not an option, check spelling or if string')
end
end