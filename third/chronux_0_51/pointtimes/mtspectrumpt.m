function [S,f,R,Serr]=mtspectrumpt(data,params,fscorr,t)
% Multi-taper spectrum - point process times
%
% Usage:
%
% [S,f,R,Serr]=mtspectrumpt(data,params,fscorr,t)
% Input: 
%       data        (structure array of spike times with dimension channels/trials; also accepts 1d array of spike times) -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over channels/trials when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
%       t        (time grid over which the tapers are to be calculated:
%                      this argument is useful when calling the spectrum
%                      calculation routine from a moving window spectrogram
%                      calculation routine). If left empty, the spike times
%                      are used to define the grid.
% Output:
%       S       (spectrum with dimensions frequency x channels/trials if trialave=0; dimension frequency if trialave=1)
%       f       (frequencies)
%       R       (rate)
%       Serr    (error bars) - only if err(1)>=1
if nargin < 1; error('Need data'); end;
if nargin < 2; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargout > 3 && err(1)==0; error('cannot compute error bars with err(1)=0; change params and run again'); end;
if nargin < 3 || isempty(fscorr); fscorr=0;end;
if nargin < 4 || isempty(t);
   [mintime,maxtime]=minmaxsptimes(data);
   dt=1/Fs; % sampling time
   t=mintime-dt:dt:maxtime+dt; % time grid for prolates
end;
N=length(t); % number of points in grid for dpss
nfft=2^(nextpow2(N)+pad); % number of points in fft of prolates
[f,findx]=getfgrid(Fs,nfft,fpass); % get frequency grid for evaluation
tapers=dpsschk(tapers,N,Fs); % check tapers
[J,Msp,Nsp]=mtfftpt(data,tapers,nfft,Fs,t,f,findx); % mt fft for point process times
S=squeeze(mean(conj(J).*J,2));
if trialave; S=squeeze(mean(S,2));Msp=mean(Msp);end;
R=Msp*Fs;
if nargout==4;
   if fscorr==1;
      Serr=specerr(S,J,err,trialave,Nsp);
   else
      Serr=specerr(S,J,err,trialave);
   end;
end;
