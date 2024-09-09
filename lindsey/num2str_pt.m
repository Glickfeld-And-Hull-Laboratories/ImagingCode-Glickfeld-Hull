function s = num2str_pt(x, p)

s = num2str(x);
ind = find(s == '.');
s = [s(1:ind-1) 'pt' s(ind+1:end)];