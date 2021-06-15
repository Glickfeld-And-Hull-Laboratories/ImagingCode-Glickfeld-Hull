function conStruct_all = nakaRushtonFit(resp_mat, cons)
    %resp_mat is nCells x nContrasts
    nCellsTot = size(resp_mat,1);
    nCon = length(cons);
    s4 = zeros(1,nCon);
    s = zeros(1);
    conStruct_all = struct('resp',s4,'fit',s4,'C50r',s,'Rsq',s,'x0',s4);
    conStruct_all(nCellsTot) = conStruct_all;
    conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
    conRng = 0:0.001:1;
    opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'

    for iCell = 1:nCellsTot
        conStruct_all(iCell).resp = resp_mat(iCell,:);
        cRi = conStruct_all(iCell).resp;
        lb = [0 0 0.1 1];
        ub = [Inf Inf 0.8 Inf];
        SStot = sum((cRi-mean(cRi)).^2);
        R2best = -Inf;
        for i=1%1:4
            x0 = [cRi(1) max(cRi) 0.1+0.1*i 3]; %BL Rmax C50 n
            [cF, res] = lsqcurvefit(conModelH,x0,cons,cRi,lb,ub,opts);
            R2 = 1-res/SStot;
            if R2>R2best
                R2best = R2;
                cFbest = cF;
                x0best = x0;
            end
        end
        cF = cFbest;
        R2 = R2best;

        conStruct_all(iCell).fit = cF;
        conStruct_all(iCell).Rsq = R2;
        conStruct_all(iCell).x0 = x0best;

        fitout = conModelH(cF,conRng);
        R50 = fitout(1)+(fitout(end)-fitout(1))/2;
        fitout50rect = abs(fitout - R50);
        i50 = find(fitout50rect == min(fitout50rect),1);
        C50 = conRng(i50);
        conStruct_all(iCell).C50r = C50;
    end