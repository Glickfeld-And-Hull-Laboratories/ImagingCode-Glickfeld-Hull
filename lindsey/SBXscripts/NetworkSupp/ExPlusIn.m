function [AmpA, AmpB,AmpC,AmpD, fastRiseTau, slowRiseTau,fastDecayTau,slowDecayTau,sse,R_square] = ExPlusIn(t,data)
    AmpA_guess = 0.5.*ones(20,1);
    AmpB_guess = 0.5.*ones(20,1);
    AmpC_guess = 0.5.*ones(20,1);
    AmpD_guess = 0.5.*ones(20,1);
    fastRiseTau_guess = 3;
    slowRiseTau_guess = 10;
    fastDecayTau_guess = 7;
    slowDecayTau_guess = 30;

    AmpA_lb = zeros(20,1);
    AmpB_lb = zeros(20,1);
    AmpC_lb = zeros(20,1);
    AmpD_lb = zeros(20,1);
    fastRiseTau_lb = 1;
    slowRiseTau_lb = 3;
    fastDecayTau_lb = 4;
    slowDecayTau_lb = 10;

    AmpA_ub = 5.*ones(20,1);
    AmpB_ub = 5.*ones(20,1);
    AmpC_ub = 5.*ones(20,1);
    AmpD_ub = 5.*ones(20,1);
    fastRiseTau_ub = 10;
    slowRiseTau_ub = 30;
    fastDecayTau_ub = 30;
    slowDecayTau_ub = 100;

    options = optimset('MaxFunEvals',inf,'MaxIter',100000);
    [out,fval,success] = fminsearchbnd(@ExPlusIn,...
        [AmpA_guess, AmpB_guess,AmpC_guess,AmpD_guess, fastRiseTau_guess, slowRiseTau_guess,fastDecayTau_guess,slowDecayTau_guess],...
        [AmpA_lb, AmpB_lb,AmpC_lb,AmpD_lb, fastRiseTau_lb, slowRiseTau_lb,fastDecayTau_lb,slowDecayTau_lb],...
        [AmpA_ub, AmpB_ub,AmpC_ub,AmpD_ub, fastRiseTau_ub, slowRiseTau_ub,fastDecayTau_ub,slowDecayTau_ub], options);

    if success == 1
        AmpA=out(1);
        AmpB=out(2);
        AmpC=out(3);
        AmpD=out(4);
        fastRiseTau=out(5);
        slowRiseTau=out(6);
        fastDecayTau=out(7);
        slowDecayTau=out(8);
        sse=fval;
        sse_tot = sum((data - mean(data)).^2);
        R_square = 1-(sse/sse_tot);
    else
        AmpA=nan;
        AmpB=nan;
        AmpC=nan;
        AmpD=nan;
        fastRiseTau=nan;
        slowRiseTau=nan;
        fastDecayTau=nan;
        slowDecayTau=nan;
        sse=nan;
        sse_tot = nan;
        R_square = nan;
    end

    function miaosse = ExPlusIn(in)
        % pull out the slope and intercept
        AmpA = in(1);
        AmpB = in(2);
        AmpC = in(3);
        AmpD = in(4);
        fastRiseTau = in(5);
        slowRiseTau = in(6);
        fastDecayTau = in(7);
        slowDecayTau = in(8);

        allBuf = zeros(1,62);
        Ibuf = zeros(1,2);
        fastI = [allBuf -exp(-t./fastDecayTau)+exp(-t./fastRiseTau)];
        fastI = [Ibuf fastI(1:end-2)];
        fastE = [allBuf exp(-t./fastDecayTau)-exp(-t./fastRiseTau)];
        Ibuf = zeros(1,4);
        slowI = [allBuf -exp(-t./slowDecayTau)+exp(-t./slowRiseTau)];
        slowI = [Ibuf slowI(1:end-4)];
        Ibuf = zeros(1,3);
        slowE = [allBuf exp(-t./slowDecayTau)-exp(-t./slowRiseTau)];
        slowE = [Ibuf slowE(1:end-3)];
        y_fit = AmpA.*fastE+AmpB.*fastI+AmpC.*slowE+AmpD.*slowI;

        residuals = data - y_fit;
        miaosse = sum(residuals.^2);
    end
end