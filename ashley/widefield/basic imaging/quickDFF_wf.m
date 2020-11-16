function data_dff = quickDFF_wf(data)

data = double(data);
data_sub = data - min(data(:));
f = mean(data_sub,3);
data_dff = (data_sub-f)./f;

end