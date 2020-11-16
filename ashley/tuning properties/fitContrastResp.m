function [fit,R2,c50] = fitContrastResp(avgResp,contrasts)
    
    conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
    conRng = 0.001:0.001:1;
    opts = optimoptions('lsqcurvefit','Display','off');
    lb = [0 0 0.1 1];
    ub = [Inf Inf 0.8 Inf];

%     nBaseFr = round(frameRateHz); % 1 second
%     nStimFr = round(frameRateHz); % 1 second
    
    [nstim,nc] = size(avgResp);
    fit = nan(length(conRng),nc);
    R2 = nan(1,nc);
    c50 = nan(1,nc);
    for icell = 1:nc
        cRi = avgResp(:,icell)';
        x0 = [cRi(1) max(cRi) 0.2 3]; %BL Rmax C50 n
        [cF, res] = lsqcurvefit(conModelH,x0,contrasts,cRi,lb,ub,opts);
        SStot = sum((cRi-mean(cRi)).^2);
        R2(:,icell) = 1-res/SStot;
        fit(:,icell) = conModelH(cF,conRng);
        R50 = fit(1,icell)+(fit(end,icell)-fit(1,icell))/2;
        i50 = find(abs(fit(:,icell) - R50) == min(abs(fit(:,icell) - R50)),1);
        c50(:,icell) = conRng(i50);
    end
end