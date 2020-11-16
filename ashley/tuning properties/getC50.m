function c50 = getC50(fit)
conRng = 0.001:0.001:1;
c50 = nan(1,size(fit,2));
for icell = 1:size(fit,2)
    R50 = fit(1,icell)+(fit(end,icell)-fit(1,icell))/2;
    i50 = find(abs(fit(:,icell) - R50) == min(abs(fit(:,icell) - R50)),1);
    c50(:,icell) = conRng(i50);
end
end