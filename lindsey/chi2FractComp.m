function [tbl ch2stat,p] = chi2FractComp(n1,N1,n2,N2);
    % n1- numerator 1
    % N1 - denominator 1
    % n2- numerator 2
    % N2 - denominator 2
    x1 = [repmat('a',N1,1); repmat('b',N2,1)];
    x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
    [tbl ch2stat,p] = crosstab(x1,x2);
end

