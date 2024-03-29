function [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = CellsortPCA2(rt, flims, nPCs, dsamp, outputdir, fname, badframes)
% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
%
% CELLSORT
% Read TIFF movie data and perform singular-value decomposition (SVD)
% dimensional reduction.
%
% Inputs:
%   fn - movie file name. Must be in TIFF format.
%   flims - 2-element vector specifying the endpoints of the range of
%   frames to be analyzed. If empty, default is to analyze all movie
%   frames.
%   nPCs - number of principal components to be returned
%   dsamp - optional downsampling factor. If scalar, specifies temporal
%   downsampling factor. If two-element vector, entries specify temporal
%   and spatial downsampling, respectively.
%   outputdir - directory in which to store output .mat files
%   badframes - optional list of indices of movie frames to be excluded
%   from analysis
%
% Outputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y array of N spatial signal mixtures sampled at
%   X x Y spatial points.
%   CovEvals - largest eigenvalues of the covariance matrix
%   covtrace - trace of covariance matrix, corresponding to the sum of all
%   eigenvalues (not just the largest few)
%   movm - average of all movie time frames at each pixel
%   movtm - average of all movie pixels at each time frame, after
%   normalizing each pixel deltaF/F
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

tic


%-----------------------
c = 'double';

global nt useframes jj jjind
useframes = setdiff((flims(1):flims(2)), badframes);
nt = length(useframes);

if nargin<3 || isempty(nPCs)
    nPCs = min(150, nt);
end
if nargin<4 || isempty(dsamp)
    dsamp = [1,1];
end
if nargin<5 || isempty(outputdir)
    outputdir = [pwd,'/cellsort_preprocessed_data/'];
end
if nargin<6
    badframes = [];
end
if isempty(dir(outputdir))
    mkdir(pwd, '/cellsort_preprocessed_data/')
end
if outputdir(end)~='/';
    outputdir = [outputdir, '/'];
end
[pixw,pixh] = size(rt(:,:,1));
% Downsampling
if length(dsamp)==1
    dsamp_time = dsamp(1);
    dsamp_space = 1;
else
    dsamp_time = dsamp(1);
    dsamp_space = dsamp(2); % Spatial downsample
%     [pixw,pixh] = size(imresize( rt(:,:,1), 1/dsamp_space, 'bilinear' ));
    [pixw,pixh] = size(imresize( rt(:,:,1), 1/dsamp_space));
end


nt = nt/dsamp_time;
npix = pixw*pixh;

fprintf('   %d pixels x %d time frames;', npix, nt)
if nt<npix
    fprintf(' using temporal covariance matrix.\n')
else
    fprintf(' using spatial covariance matrix.\n')
end

% Create covariance matrix
if nt < npix
    [covmat, mov, movm, movtm] = create_tcov(rt, pixw, pixh, useframes, nt, dsamp, c);
else
    [covmat, mov, movm, movtm] = create_xcov(rt, pixw, pixh, useframes, nt, dsamp, c);
end
covtrace = trace(covmat) / npix;
movm = reshape(movm, pixw, pixh);

if nt < npix
    % Perform SVD on temporal covariance
    [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix, c);

    % Load the other set of principal components
    [mixedfilters] = reload_moviedata(pixw*pixh, mov, mixedsig, CovEvals);
else
    % Perform SVD on spatial components
    [mixedfilters, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix, c);

    % Load the other set of principal components
    [mixedsig] = reload_moviedata(nt, mov', mixedfilters, CovEvals);
end
mixedfilters = reshape(mixedfilters, pixw,pixh,nPCs);
fnmat = [outputdir, fname, '_','PCA', '_', date,'.mat'];

firstframe_full = rt(:,:,1);
firstframe = firstframe_full;
if dsamp_space>1
    firstframe = imresize(firstframe, size(mov(:,:,1)),'bilinear');
end

%------------

save(fnmat,'mixedfilters','CovEvals','mixedsig', ...
    'movm','movtm','covtrace')
fprintf(' CellsortPCA: saving data and exiting; ')

toc

    function [covmat, mov, movm, movtm] = create_xcov(rt, pixw, pixh, useframes, nt, dsamp, sclass)
        %-----------------------
        % Load movie data to compute the spatial covariance matrix

        npix = pixw*pixh;

        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

        if (dsamp_space==1)
            mov = zeros(pixw, pixh, nt, sclass);
            for jjind=1:nt
                jj = useframes(jjind);
                mov(:,:,jjind) = rt(:,:,jj);
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        else
%             [pixw_dsamp,pixh_dsamp] = size(imresize( rt(:,:,1), 1/dsamp_space, 'bilinear' ));
            [pixw_dsamp,pixh_dsamp] = size(imresize( rt(:,:,1), 1/dsamp_space));
            mov = zeros(pixw_dsamp, pixh_dsamp, nt, sclass);
           
            for jjind=1:nt
                jj = useframes(jjind);
%                 mov(:,:,jjind) = imresize( rt(:,:,jj), 1/dsamp_space, 'bilinear' );
                mov(:,:,jjind) = imresize( rt(:,:,jj), 1/dsamp_space);
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movm = max(movm,1);
        movmzero = (movm==0);
        movm(movmzero) = 1;

        mov = mov ./ (movm * ones(1,nt)) - 1; % Compute Delta F/F
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp_time,1)/dsamp_time, 1, mov, [], 2);
            mov = downsample(mov', dsamp_time)';
        end

        movtm = mean(mov,2); % Average over space
        clear movmzeros

        c1 = (mov*mov')/size(mov,2);
        toc
        covmat = c1 - movtm*movtm';
        clear c1
    end

    function [covmat, mov, movm, movtm] = create_tcov(rt, pixw, pixh, useframes, nt, dsamp, sclass)
        %-----------------------
        % Load movie data to compute the temporal covariance matrix
        npix = pixw*pixh;

        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

        if (dsamp_space==1)
            mov = zeros(pixw, pixh, nt, sclass);
            for jjind=1:nt
                jj = useframes(jjind);
                mov(:,:,jjind) = rt(:,:,jj);
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        else
%             [pixw_dsamp,pixh_dsamp] = size(imresize( rt(:,:,1), 1/dsamp_space, 'bilinear' ));
            [pixw_dsamp,pixh_dsamp] = size(imresize( rt(:,:,1), 1/dsamp_space));
            mov = zeros(pixw_dsamp, pixh_dsamp, nt, sclass);
            
            
            for jjind=1:nt
                jj = useframes(jjind);
%                 mov(:,:,jjind) = imresize( rt(:,:,jj), 1/dsamp_space, 'bilinear' );
                mov(:,:,jjind) = imresize( rt(:,:,jj), 1/dsamp_space);
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0); % Avoid dividing by zero
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1;
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp_time,1)/dsamp_time, 1, mov, [], 2);
            mov = downsample(mov', dsamp_time)';
        end

        c1 = (mov'*mov)/npix;
        movtm = mean(mov,1); % Average over space
        covmat = c1 - movtm'*movtm;
        clear c1
    end

    function [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix, sclass)
        %-----------------------
        % Perform SVD

        covtrace = trace(covmat) / npix;

        opts.disp = 0;
        opts.issym = 'true';
        if nPCs<size(covmat,1)
            if strcmp(sclass, 'gpuArray')
                covmat = gather(covmat);
                [mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
                mixedsig = gpuArray(mixedsig);
                CovEvals = gpuArray(CovEvals);
            else
                [mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
            end
         
        else
            [mixedsig, CovEvals] = eig(covmat);
            CovEvals = diag( sort(diag(CovEvals), 'descend'));
            nPCs = size(CovEvals,1);
        end
        CovEvals = diag(CovEvals);
        if nnz(CovEvals<=0)
            nPCs = nPCs - nnz(CovEvals<=0);
            fprintf(['Throwing out ',num2str(nnz(CovEvals<0)),' negative eigenvalues; new # of PCs = ',num2str(nPCs),'. \n']);
            mixedsig = mixedsig(:,CovEvals>0);
            CovEvals = CovEvals(CovEvals>0);
        end

        mixedsig = mixedsig' * nt;
        CovEvals = CovEvals / npix;

        percentvar = 100*sum(CovEvals)/covtrace;
        fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
    end

    function [mixedfilters] = reload_moviedata(npix, mov, mixedsig, CovEvals)
        %-----------------------
        % Re-load movie data
        nPCs = size(mixedsig,1);

        Sinv = inv(diag(CovEvals.^(1/2)));

        movtm = mean(mov,1); % Average over space
        movuse = mov - ones(npix,1) * movtm;
        mixedfilters = reshape(movuse * mixedsig' * Sinv, npix, nPCs);
    end

    function j = tiff_frames(rt)
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
                rt(:,:,j);
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
    end
end