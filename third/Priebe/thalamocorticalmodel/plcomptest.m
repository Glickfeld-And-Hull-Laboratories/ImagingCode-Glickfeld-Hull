function [Rp,Rc] = plcomptest(gratedirs,plaiddirs,theta)
%
%  function [Rp,Rc] = plcomptest(gratedirs,plaiddirs,theta)
%
%


if ~(exist('theta'))
  theta = 135;
end

ndirs = length(gratedirs);

nrots = (theta * ndirs/360)/2;

comp_pred = circshift(gratedirs,nrots) + circshift(gratedirs,-nrots);
pat_pred = gratedirs;
c = corrcoef(plaiddirs,pat_pred);
corrp = c(1,2);
c = corrcoef(plaiddirs,comp_pred);
corrc = c(1,2);
c = corrcoef(comp_pred,pat_pred);
corrpred = c(1,2);

Rp = (corrp - corrc*corrpred) / (sqrt((1-corrc^2)*(1-corrpred^2)));

Rc = (corrc - corrp*corrpred) / (sqrt((1-corrp^2)*(1-corrpred^2)));
