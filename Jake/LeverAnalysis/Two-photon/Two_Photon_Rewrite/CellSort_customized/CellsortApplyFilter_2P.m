function cell_sig = CellsortApplyFilter_2P(fn, ica_segments, flims, movm, subtractmean, normal)
% cell_sig = CellsortApplyFilter(fn, ica_segments, flims, movm, subtractmean)
%
%CellsortApplyFilter
% Read in movie data and output signals corresponding to specified spatial
% filters
%
% Inputs:
%     fn - file name of TIFF movie file
%     ica_segments - nIC x X matrix of ICA spatial filters
%     flims - optional two-element vector of frame limits to be read
%     movm - mean fluorescence image
%     subtractmean - boolean specifying whether or not to subtract the mean
%     fluorescence of each time frame
%
% Outputs:
%     cell_sig - cellular signals
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

if (nargin<3)||isempty(flims)
    nt = size(fn,3);
    flims = [1,nt];
else
    nt = diff(flims)+1;
end
if nargin<5
    subtractmean = 1;
end

[pixw,pixh,~] = size(fn);
if (nargin<4)||isempty(movm)
    movm = ones(pixw,pixh);
else
    movm = double(movm);
end
movm(movm==0) = 1; % Just in case there are black areas in the average image
k=0;

cell_sig = zeros(size(ica_segments,1), nt);
ica_segments = reshape(ica_segments, [], pixw*pixh);

fprintf('Loading %5g frames for %d ROIs.\n', nt, size(ica_segments,1))
while k<nt
    ntcurr = min(500, nt-k);
    mov = zeros(pixw, pixh, ntcurr);
    for j=1:ntcurr
        movcurr = fn(:,:, j+k+flims(1)-1);
        mov(:,:,j) = movcurr;
    end
    if normal == 1
        mov = (mov ./ repmat(movm, [1,1,ntcurr])) - 1; % Normalize by background and subtract mean
    end
    
    if subtractmean
        % Subtract the mean of each frame
        movtm = mean(mean(mov,1),2);
        mov = mov - repmat(movtm,[pixw,pixh,1]);
    end

    mov = reshape(mov, pixw*pixh, ntcurr);
    cell_sig(:, k+[1:ntcurr]) = ica_segments*mov;

    k=k+ntcurr;
    fprintf('Loaded %3.0f frames; ', k)
    toc
end

function j = tiff_frames(fn)
%
% n = tiff_frames(filename)
%
% Returns the number of slices in a TIFF stack.
%
%

status = 1; j=0;
jstep = 10^3;
while status
    try
        j=j+jstep;
        imread(fn,j);
    catch
        if jstep>1
            j=j-jstep;
            jstep = jstep/10;
        else
            j=j-1;
            status = 0;
        end
    end
end