function [n, n2] = subplotn(ntot);
%output minimum subplot dimensions given n conditions
    if ntot<4
        n = ntot;
        n2 = 1;
    else
        n = ceil(sqrt(ntot));
        if n^2-n >= ntot
            n2 = n-1;
        else
            n2 =n;
        end
    end
end