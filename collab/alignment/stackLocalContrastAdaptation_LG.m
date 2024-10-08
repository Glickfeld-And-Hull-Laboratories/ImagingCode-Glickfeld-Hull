function out = stackLocalContrastAdaptation(stack,sigma,offset)
%STACKLOCALCONTRASTADAPTATION 
% OUT = STACKLOCALCONTRASTADAPTATION(STACK,SIGMA,OFFSET) where SIGMA is the
% size of the smoothing kernel (in pixels) and OFFSET is the size of the
% offset added to avoid dividing by zero. OFFSET is expressed in percent of
% the range of the input stack.
%
% by Vincent Bonin 
% based on stack_localcontrastadj by Mark histed

% Change log
% 03/18/08 VB Started from stack_localcontrastadj
%             Streamlined, removed type cast. Additional arguments can be
%             put back in as needed.
% 03/23/08 VB Tested vs. stack_localcontrastadj and very similar images.
%

[ny,nx,nframes]=size(stack);

av = mean(stack,3);

% lowspass filter
h = fspecial('gaussian',3*sigma*[1 1],sigma);
for n = 1:nframes;
   smFr = imfilter(stack(:,:,n),h);
   smOffset = range(smFr(:)) * offset/100;   
   smFr = smFr + smOffset;
   out=repmat(single(smFr),[1 1 nframes]);
   out(:,:,n)=stack(:,:,n)./out;
end
% add offset (to avoid divide by zero)

%{
out = zeros(ny,nx,nframes,'single');

for index = 1:nframes
    out(:,:,index) = stack(:,:,index)./smFr;
end
%}









return;
