function out = stackGroupProject_timecourses(array,ratio)
%STACKGROUPPROJECT Downsample stack array by ratio
%   out = stackGroupProject(array,ratio)

[nframes,nx]=size(array);

% which frames?
grpfrns = floor(nframes/ratio)*ratio; 
if grpfrns ~= nframes
    % truncate
    desfrns = 1:grpfrns;
else
    desfrns = 1:nframes;
end
ndes = length(desfrns);

siz = ndes / ratio;
%{
c = class(array);
out = zeros(ny,nx,siz,c);

for i = 1:siz
    out(:,:,i)=mean(array(:,:,(i-1)*ratio+1:i*ratio),3);
end
%}

array=squeeze(mean(reshape(array(desfrns,:),ratio,siz,nx),1));

out=array;

return