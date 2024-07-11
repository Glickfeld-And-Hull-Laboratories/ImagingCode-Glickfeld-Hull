function [outs,stack]=stackRegister(stack,target,usFac);
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
% [OUTS]=STACKREGISTER(STACK,TARGET,USFAC) to oversample (default 100)
%
% [OUTS,REG]=STACKREGISTER(STACK,TARGET) returns registered stack REG.
%
% based on DFTREGISTRATION by Manuel Guizar
%
% Notes: For MATLAB 2008b and earlier:
%          In File->Preferences->Multithreading, 
%          you should set it to Manual, *1* thread.
%        For MATLAB 2009a and later, can expect 2x speedup on multicore machin
%
% todo: - function that lets user specify shifts
%       - argument that lets user specify roi over which correlation computed

% vb 09/04/17 now processes frames as singles because faster
% vb 09/05/25 sets max num threads to 2 because maginal speed up for > 2


% $$$ if ~verLessThan('matlab','7.5')
% $$$     prevn = maxNumCompThreads(1); % marginal gain beyond 1 thread
% $$$ end

if nargin < 3, usFac = 100; end
c = class(stack);

TARGET = fft2(cast(target,'single'));

[ny,nx,nframes]=size(stack);

outs = zeros(nframes,4);

% if nargout > 1
%     reg = zeros(size(stack),c);
% end

verbose_threshold = 64;

if nframes > verbose_threshold
    tic;
    fprintf(1, 'Starting, frame 0\n');
end

for index = 1:nframes
    if mod(index,250)==0 
        fprintf(1,'Frame %i (%2.1f fps)\n',index,index/toc);
    end
    SLICE = fft2(cast(stack(:,:,index),'single'));
    [outs(index,:) temp ] = dftregistration(TARGET,SLICE,usFac);

    if nargout > 1
        wS = warning('off'); 
        if strcmp(c, 'gpuArray')
            stack(:,:,index) = abs(ifft2(temp));
        else
            stack(:,:,index) = cast(abs(ifft2(temp)),c);
        end
        
        warning(wS);
    end
end

if nframes > verbose_threshold
    t = toc;
    fprintf('Registered %i frames in %2.1f seconds (%2.1f fps)\n',nframes,t,nframes/t);
end

% $$$ if ~verLessThan('matlab','7.5')    
% $$$     maxNumCompThreads(prevn);  % restore default
% $$$ end

return;
