function result = stackFilter(stack,kernel);
%STACKFILTER Filter each slice of stack
% RESULT = STACKFILTER(STACK);
% RESULT = STACKFILTER(STACK,SIGMA); 
% RESULT = STACKFILTER(STACK,KERNEL);
%
% see MATLAB's IMFILTER
% see core's SMOOTH2, IMFILTER2
%
% this function is needed becasue SMOOTH2 does not operate on uint16
%

if nargin < 2
    sigma = 3;
    kernel = fspecial('gauss',[1 1]*sigma*6,sigma);
end

if nargin > 1 & isscalar(kernel)
    sigma = kernel;
    kernel = fspecial('gauss',[1 1]*sigma*6,sigma);
end
c = class(stack);

[ny,nx,nslices]=size(stack);

result = zeros(ny,nx,nslices,c);

for islice = 1:nslices
    if mod(islice,100)==1
        fprintf(1,'slice %i\n',islice);
    end
    result(:,:,islice)=imfilter(stack(:,:,islice),kernel,'conv');
end
    
return;

