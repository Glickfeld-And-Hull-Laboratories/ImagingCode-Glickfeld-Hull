%function eegpower


%Fs = 500; % Sampling frequency
%T = 1/Fs; % Sample time
%L = 4000; % Length of signal
%t = (0:L-1)*T; % Time vector

%NFFT = 2^nextpow2(L); % Next power of 2 from length of y X = fft(x,NFFT)/L;
%f = Fs/2*linspace(0,1,NFFT/2+1);
%AMP = 2*abs(X(1:NFFT/2+1));

%Fs = 500; % Sampling frequency 
%T = 1/Fs; % Sample time 
%L = 4000; % Length of signal 
%t = (0:L-1)*T; % Time vector
%x = cos(2*pi*100*t)+randn(size(t));
%xdft = fft(x);
%Pxx = 1/(L*Fs)*abs(xdft(1:length(x)/2+1)).^2;
%freq = 0:Fs/L:Fs/2;
%plot(freq,10*log10(Pxx));
%xlabel('Hz'); ylabel('dB/Hz');
%Pxx(2:end-1) = 2*Pxx(2:end-1);
%plot(freq,10*log10(Pxx));


data=double(NS5.Data(1,:));
figure(1); 
xlabel('sample'); 
ylabel('magnitude'); 
plot(data); 
legend('raw eeg plot'); 
%data= data'; 
%value1=data; 
chan1= data-mean(data); 
fs=30000;
d=1/fs; 
t= 0:0.0001:1-0.0001 %t=[0:length(data)-1]*d; 
figure(2); 
plot(chan1);
title('original signal'); 
fs=fft(chan1,10000); 
pp=abs(fs)/128; %pp=fs.*conj(fs)/10000; 
L=length(fs)/2;
ff=(0:L)/10000/d; 
figure(3); 
plot(ff,pp(0:4999)); 
ylabel('power spectrum density');
xlabel('frequency');
title('signal power spectrum');