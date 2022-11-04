function x = gengauss(w,p,c)
%GENGAUSS Generate gaussian peak of unity height.
%   X = GENGAUSS(W)  or   X = GENGAUSS(W,P,C)
%
%    where:
%
%    W    is the distance from center to 2% of full height
%    P    (optional) is the number of points x will contain
%    C    (optional) is the point number for the peak center
%
%    The peak will be centered in the vector unless P and C are specified.
%
%

%   by Richard Kramer.
%   Copyright (c) 1988-1993 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 1993/09/01 19:37:08 $

if nargin == 2, c = 0; end
if nargin == 1, c = 0; p = 0; end

if  p < 1, points = round(3 * w);
    else points = round(p);
end

if c < 1, center = round(points/2);
    else center = round(c);
end

w = w/2;
if c < 0, center = center - c; end

for n = 1:points,
    x(n) = exp(-((n-center)/w)^2);
end

