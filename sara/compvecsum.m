function [amp,ang] = compvecave(amps, angs)
%  function [amp,ang] = compvecave(amps, angs)
comprep = zeros(size(amps,1),1);

for j = 1:size(amps,1)
  comprep(j) = cos(angs(j))*amps(j) + i*sin(angs(j))*amps(j);
end
amp = abs(sum(comprep));
ang = angle(sum(comprep));
