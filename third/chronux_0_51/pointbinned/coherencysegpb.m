function [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpb(data1,data2,win,params,segave,fscorr)
% Multi-taper coherency,cross-spectrum and individual spectra computed by segmenting two univariate binned point processes into chunks 
%
% Usage:
% [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpb(data1,data2,win,params,segave,fscorr)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data1 (column vector, binned point process data) -- required
%       data2 (column vector, binned point process data) -- required
%       win   (length of segments) - required
%       params: structure with fields tapers, pad, Fs, fpass, err
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
%       segave (average over segments for 1, don't average for 0)
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       C (magnitude of coherency - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       phi (phase of coherency - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S12 (cross spectrum -  frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S1 (spectrum 1 - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S2 (spectrum 2 - frequencies x segments if segave=0; dimension
%       frequencies if segave=1)
%       f (frequencies)
%       zerosp (1 for segments where no spikes were found, 0 otherwise)
%       confC (confidence level for C at 1-p %)
%       phierr (error bars for phi)
%       Cerr  (Jackknife error bars for C - use only for Jackknife)

if nargin < 3; error('Need data1 and data2 and size of segment'); end;
if nargin < 4; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 5 || isempty(segave); segave=1;end;
if nargin < 6 || isempty(fscorr); fscorr=0; end;

if nargout > 9 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
if nargout > 7 && err(1)==0;
    error('When error are desired, err(1) has to be non-zero.');
end;

[N1,C1,N2,C2]=check_consistency(data1,data2);
N=N1;
dt=1/Fs; % sampling interval
T=N*dt; % length of data in seconds
E=0:win:T-win; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
data1=createdatamatpb(data1,E,Fs,win); % segmented data 1
data2=createdatamatpb(data2,E,Fs,win); % segmented data 2
params.trialave=segave;
if err(1)==0;
   [C,phi,S12,S1,S2,f,zerosp]=coherencypb(data1,data2,params,fscorr); % compute coherency for segmented data
elseif err(1)==1;
   [C,phi,S12,S1,S2,f,zerosp,confC,phierr]=coherencypb(data1,data2,params,fscorr); % compute coherency for segmented data
elseif err(1)==2;
   [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypb(data1,data2,params,fscorr); % compute coherency for segmented data
end;
