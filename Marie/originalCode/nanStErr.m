function sterr = nanStErr(A)
stdev = nanstd(A);
n = length(isfinite(A));
sterr = stdev/sqrt(n);
end
