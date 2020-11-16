function [lowerErr, upperErr, lowerCI, upperCI] = ciFromBoot(data,ciPct)
%%% bootstraps should be in first dimension of data
%%% ciPct should be between 1 to 99;

nBoot = size(data,1);
pctEndRange = round((100-ciPct)/2/100*nBoot);

if pctEndRange == 0
    ind_lower = 1;
else
    ind_lower = pctEndRange;
end
ind_upper = nBoot-pctEndRange;

data_sort = sort(data,1);
data_mean = mean(data,1);

lowerCI = data_sort(ind_lower,:);
upperCI = data_sort(ind_upper,:);

lowerErr = data_mean - lowerCI;
upperErr = upperCI - data_mean;

end