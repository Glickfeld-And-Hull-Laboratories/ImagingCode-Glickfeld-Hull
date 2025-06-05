function time = make_time(signal,Fs,dim)

time = 0:1/Fs:(1/Fs*size(signal,dim))-1/Fs;
end