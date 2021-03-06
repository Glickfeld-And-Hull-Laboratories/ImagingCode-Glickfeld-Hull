function [output] = usealgebraicSB(model)
% usealgebraicSB: checks if an SBmodel contains algebraic rules.
%
% Output Arguments:
% =================
% output: =1 if algebraic rules are used, =0 if not

% Information:
% ============
% Copyright (C) 2008  Henning Schmidt, henning@sbtoolbox2.org
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% USA.

if ~isSBmodel(model),
    error('Given model is not an SBmodel.');
end

% default setting
output = 0; % no algebraic rule present

% get model structure
ms = struct(model);

% check algebraic rules
if ~isempty(ms.algebraic),
    output = 1;
end
