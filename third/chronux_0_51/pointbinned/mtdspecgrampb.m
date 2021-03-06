function [dS,t,f]=mtdspecgrampb(data,movingwin,phi,params)
% Multi-taper derivatives of time-frequency spectrum - binned point process
%
% Usage:
%
% [dS,t,f]=mtdspecgrampb(data,movingwin,phi,params)
% Input: 
%   Note that all times can be in arbitrary units. But the units have to be
%   consistent. So, if E is in secs, win, t have to be in secs, and Fs has to
%   be Hz. If E is in samples, so are win and t, and Fs=1. In case of spike
%   times, the units have to be consistent with the units of data as well.
%       data        (in form samples x channels/trials) -- required
%       movingwin         (in the form [window winstep] i.e length of moving
%                                                 window and step size.
%                                                 Note that units here have
%                                                 to be consistent with
%                                                 units of Fs
%       phi         (angle for evaluation of derivative) -- required.
%                       e.g. phi=[0,pi/2] giving the time and frequency
%                       derivatives
%       params: structure with fields tapers, pad, Fs, fpass,trialave
%       -optional
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and
%                                   Fs/2
%           trialave (average over trials when 1, don't average when 0) -
%           optional. Default 0
% Output:
%       dS      (spectral derivative in form phi x time x frequency x channels/trials if trialave=0; phi x time x frequency if trialave=1)
%       t       (times)
%       f       (frequencies)

if nargin < 3; error('Need data, window parameters and angle'); end;
if nargin < 4; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);


N=size(data,1);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=2^(nextpow2(Nwin)+pad);
[f,findx]=getfgrid(Fs,nfft,fpass); 
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers

winstart=1:Nstep:N-Nwin+1;
nw=length(winstart);
for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   [ds,f]=mtdspectrumpb(datawin,phi,params);
   dS(n,:,:,:)=ds;
end;
sz=size(ds);
if length(sz)==3;
   dS=permute(dS,[2 1 3 4]);
else
   dS=permute(dS,[2 1 3]);
end;
winmid=winstart+round(Nwin/2);
t=winmid/Fs;
