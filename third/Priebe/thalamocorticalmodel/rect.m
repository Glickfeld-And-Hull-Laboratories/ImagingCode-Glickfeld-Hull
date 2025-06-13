function mm = rect(inmm)
mm = zeros(size(inmm));
mm(find(inmm>0)) = inmm(find(inmm>0));
