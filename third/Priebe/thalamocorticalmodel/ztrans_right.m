function [vals ] = ztrans(invals,ndir)
vals = 0.5*log( (1+invals(:,1))./ (1-invals(:,1)))/sqrt(1/(ndir-3));;
