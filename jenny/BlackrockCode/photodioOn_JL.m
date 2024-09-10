function [stimOn,stimOff,waitflag] = photodioOn_JL(Triggers,Fs,thresholdFactor,dio_idx,varargin)

p = inputParser;
p.addParamValue('stimTime', 0.2, @isnumeric);

% parse inputs
p.parse(varargin{:});
params = p.Results;

if isempty(dio_idx)
dio_idx = Triggers.photoID;
else
    waitflag = dio_idx;
end
photo_dio = Triggers(dio_idx,:);


fpass  =    [80 400];  % pass freq
fstop  =    [70 500];  % stop freq outside pass
Rpass  =    0.5;       % Attenuation(dB) for passband
Astop  =    30;        % Attenuation(dB) for stopband
n      =    cheb2ord(fpass/Fs*2,fstop/Fs*2,Rpass,Astop);   % order of chey filter

[z,p,k] =   cheby2(n,Astop,fstop/Fs*2);   % zeros, poles, and gain
[s,g]   =   zp2sos(z,p,k);                  % create second order section
Hd      =   dfilt.df2sos(s,g);                 % dfilt object


photo_dio =    filtfilthd(Hd,photo_dio);    % apply filter

photo_dio = bsxfun(@rdivide,photo_dio,(max(photo_dio,[],2)));
is_dio_on = photo_dio(500:end) > thresholdFactor * (max(photo_dio(:)));

temp_stim = find(is_dio_on) + 500;
stim_tDiff = params.stimTime;
temp_stim2 = diff([-10000 temp_stim])>stim_tDiff*Fs;
temp_stim3 = logical([diff(temp_stim)>stim_tDiff*Fs 1]);

stimOn = temp_stim(temp_stim2);
stimOff = temp_stim(temp_stim3);

if numel(stimOn) > numel(stimOff)
    stimOn = stimOn(1:numel(stimOff));
end
