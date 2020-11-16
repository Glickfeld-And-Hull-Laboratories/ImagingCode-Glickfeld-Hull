function [data,nframes] = loadsbx_choosepmt(pmt,mouse,expdate,imgfolder,fName,varargin)
%% varargin has two spots for variables, 1) n frames 2) other data path
%% Set current directory to imaging data location
if length(varargin) > 1
%     rc = varargin{2};
%     fdir = fullfile(rc.ashleyData,mouse,'two-photon imaging',expdate,imgfolder);
    fdir = varargin{2};
else
    fdir = fullfile('Z:\home\ashley\data', mouse, 'two-photon imaging',expdate,imgfolder);
end
cd(fdir);

%% load sbx info file
imgMatFile = [fName '.mat'];
load(imgMatFile);

%% get number of frames
if ~isempty(varargin)    
    if isnan(varargin{1})
        nframes = info.config.frames;
    elseif ~isempty(varargin{1})
        nframes = varargin{1}; 
    else
        nframes = info.config.frames;
    end
else
    nframes = info.config.frames;
end

%% load data
data = sbxread(fName,0,nframes);
if pmt ~= 3
    data = squeeze(data(pmt,:,:,:));
end
data = squeeze(data);
end