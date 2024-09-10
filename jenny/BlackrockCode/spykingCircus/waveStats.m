function [pt_ms,pt_ratio,pt_distance,asym,hfw,endslope]=waveStats(Wave)

Wave = interp(Wave,4); % try to get a smoother curve



[tro,troloc] = min(Wave);
[pk2,pk2loc] = max(Wave(troloc:end));
[pk1,pk1loc] = max(Wave(1:troloc));
pt_ms = (1/4)*(1/30)*pk2loc; % readjust due to the interplation for 1/4
pt_ratio = abs(pk2/tro);
pt_distance = pk2-tro; % the total amplitude in uV
asym = (pk2-pk1)/(pk2+pk1);% more close to 1 means more asymmetry, more close to 0 means more symmetrical. 
halfheight = 0.5*tro;
hfw = (1/4)*(1/30)*(find(Wave(troloc:end)>halfheight,1)+troloc-find(Wave(1:troloc)<halfheight,1)); % in miliseconds
% auc = sum((Wave(troloc:end)>0).*Wave(troloc:end))/tro;
endslope = diff(Wave([(troloc+(0.45*4*30)),(troloc+(0.55*4*30))]));
end
