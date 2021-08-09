function skipMissingPlot (yVals, ErrVals, varargin)
length(varargin)
if nargin == 2
xvals = 1:length(yVals);
xvals = xvals.';
ValidIndex = (~isnan(yVals));
xvals = xvals(ValidIndex);
ErrVals = ErrVals(ValidIndex);
yVals = yVals(ValidIndex);
xyE = [xvals, yVals, ErrVals];
errorbar(xyE(:,1), xyE(:,2), xyE(:,3));
end

if nargin == 1
xvals = 1:length(yVals);
xvals = xvals.';
ValidIndex = (~isnan(yVals));
xvals = xvals(ValidIndex);
yVals = yVals(ValidIndex);
xyE = [xvals, yVals];
plot(xyE(:,1), xyE(:,2));
end
end


