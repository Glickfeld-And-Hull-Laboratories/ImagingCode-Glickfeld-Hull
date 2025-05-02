function [outs,stack]=stackRegGPU(stack,target,usFacs,shifts_xy)
%STACKREGISTER Fourier-domain subpixel 2d rigid body registration.
% [OUTS]=STACKREGISTER(STACK,TARGET) where 
%      stack is 3d array containing stack to register
%      TARGET is 2d array 
%      OUTS(:,1) is correlation coefficients
%      OUTS(:,2) is global phase difference between images 
%               (should be zero if images real and non-negative).
%      OUTS(:,3) is net row shift
%      OUTS(:,4) is net column shift
%
% [OUTS,REG]=STACKREGISTER(STACK,TARGET) returns registered stack REG.
%
% based on DFTREGISTRATION by Manuel Guizar
% TH ADDED APR 2025: Based on stackRegister_MA, shifting operations
% chunk-wise onto the GPU. VRAM usage during is 12G max. CPU ~11% usage
% with GPU use, instead of CPU ~24% with CPU only.
% For a dataset of 529x796x86400:
% ~4.05x performance improvement for registration over stackRegister. 
% (5349s vs 1320s) 90min -> 22min
% ~2.8x performance improvement for RE-registration over stackRegister_MA
% (1317s vs 470s) 22min -> 8min
%
%MA ADD DEC 09: shifts_xy is a nframes x 4 vector identical to output of
%'outs', above, and used to corregister e.g. the red channel using green
%correg estimates
%
% Note: In File->Preferences->Multithreading, you should set it to Manual,
%    4 threads (on zquad)
%    see test results below: going from 4 to 8 threads doesn't speed it up
%    but uses more CPU
%
% todo: - function that lets user specify shifts
%       - argument that lets user specify roi over which correlation computed
% testing:
% cpu is total cpu time, summed over 8 cores, on zquad
% stack: 100 frames, 512x512, uint8
% 4 threads: cpu 28%, Registered 100 frames in 55.6 seconds (1.8)
% 6 threads: cpu 34%, Registered 100 frames in 55.9 seconds (1.8)
% 8 threads: cpu 42%, Registered 100 frames in 55.8 seconds (1.8)
%NOTE: SAME AS VB's stackregister, but option to give shifts as input
%variable

if nargin < 4, shifts_xy = []; end

if nargin < 3, usFacs = 100; end

chunkSize = 5000;
c = class(stack);

shifts_xy = gpuArray(shifts_xy);
TARGET = gpuArray(fft2(double(target)));

[ny,nx,nframes]=size(stack);

% outs = zeros(nframes,4,"gpuArray");
outs = zeros(nframes,4);

% if nargout > 1
%     reg = zeros(size(stack),c);
% end

% first registration

if ~isempty(shifts_xy)
    %save time by doing this only once
    SLICE = fft2(double(stack(:,:,1)));
    [nr,nc]=size(SLICE);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    outs = gather(shifts_xy);
end

nChunks = floor(nframes/chunkSize) + 1;
rmd = mod(nframes,chunkSize);

tic;
fprintf(1, 'Starting Registration of %i frames on the GPU\nSplit into %i chunks of %i frames with a remainder %i frames\n',nframes,nChunks-1,chunkSize,rmd);
for chk = 1:nChunks
    cStartFrame = chunkSize * (chk-1) +1;
    if chk == max(nChunks)
        cEndFrame = nframes;
        chunkNFrames = rmd;
    else
        cEndFrame = chunkSize * chk;
        chunkNFrames = chunkSize;
    end
    this_gpu_stack = gpuArray(stack(:,:,cStartFrame:cEndFrame));
    if ~isempty(shifts_xy)
        this_shifts_xy = shifts_xy(cStartFrame:cEndFrame,:);
    end
    for index = 1:chunkNFrames
        if mod(index,250)==0 
            fprintf(1,'Chunk %i/%i Frame %i/%i (%2.1f fps)\n',chk,nChunks,index,chunkNFrames,((chk-1)*chunkSize+index)/toc);
        end
        SLICE = fft2(double(this_gpu_stack(:,:,index)));
    
        %note: if shifts_xy = [], then standard correg, otherwise, uses this
        %program to return shifted version of fft2...
        if isempty(shifts_xy)
            [outs((chk-1)*chunkSize+index,:),temp] = dftregistration(TARGET,SLICE,usFacs);
        else
            row_shift = this_shifts_xy(index,3);
            col_shift = this_shifts_xy(index,4);
            diffphase = this_shifts_xy(index,2);
    
            temp = SLICE.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
            temp = temp*exp(i*diffphase);
        end
        if nargout > 1
            wS = warning('off'); 
            if class(temp) == "gpuArray"
                stack(:,:,(chk-1)*chunkSize+index) = gather(abs(ifft2(temp)));  
            else
                stack(:,:,(chk-1)*chunkSize+index) = abs(ifft2(temp));
            end
            warning(wS);
        end
    end
end
t = toc;
% outs = gather(outs);

fprintf('Registered %i frames in %2.1f seconds (%2.1f fps)\n',nframes,t,nframes/t);
gpuDevice([]);
return;
