function y = isSBML_ModifierSpeciesReference(varargin)
% isSBML_ModifierSpeciesReference(SBMLStructure, Level, Version(optional)) 
% checks that SBMLStructure represents a modifier species reference 
% within an sbml model of the specified level
% 
% if SBMLStructure represents a modifier species reference within an SBML model
% it has the appropriate fields (ONLY IN LEVEL 2)
% eg    Typecode (L2V1)
%       Metaid (L2V1)
%       Notes (L2V1)
%       Annotations (L2V1)
%       SBOTerm (L2V2)
%       Species (L2V1)
%       Id (L2V2)
%       Name (L2V2)
%
% NOTE: content of brackets indicates the level and version of sbml from which the given field
% is appropriate.
%
% Returns 1 if SBMLStructure is a structure containing each of the above
% fields (appropriate with the given level and version) 
% and the typecode is "SBML_MODIFIER_SPECIES_REFERENCE"
% 
% Returns 0 if SBMLStructure is not a structure 
% or does not contain one of the appropriate fields
% or the typecode is not "SBML_MODIFIER_SPECIES_REFERENCE"

%  Filename    :   isSBML_ModifierSpeciesReference.m
%  Description :
%  Author(s)   :   SBML Development Group <sbml-team@caltech.edu>
%  $Id: isSBML_ModifierSpeciesReference.m,v 1.9 2008/01/18 21:16:53 sarahkeating Exp $
%  $Source v $
%
%<!---------------------------------------------------------------------------
% This file is part of SBMLToolbox.  Please visit http://sbml.org for more
% information about SBML, and the latest version of SBMLToolbox.
%
% Copyright 2005-2007 California Institute of Technology.
% Copyright 2002-2005 California Institute of Technology and
%                     Japan Science and Technology Corporation.
% 
% This library is free software; you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation.  A copy of the license agreement is provided
% in the file named "LICENSE.txt" included with this software distribution.
% and also available online as http://sbml.org/software/sbmltoolbox/license.html
%----------------------------------------------------------------------- -->


%input arguments
if (nargin < 2 || nargin > 3)
    error('wrong number of input arguments');
end;

SBMLStructure = varargin{1};
Level = varargin{2};

if (nargin == 3)
    Version = varargin{3};
else
    Version = 1;
end;

if (Level == 1)
    y = 0;
    return;
else
    if (Version == 1)
        SBMLfieldnames = {'typecode', 'metaid', 'notes', 'annotation','species'};
        nNumberFields = 5;
    elseif (Version == 2)
        SBMLfieldnames = {'typecode', 'metaid', 'notes', 'annotation','species', 'id', 'name', 'sboTerm'};
        nNumberFields = 8;
    elseif (Version == 3)
        SBMLfieldnames = {'typecode', 'metaid', 'notes', 'annotation', 'sboTerm', 'species', 'id', 'name'};
        nNumberFields = 8;
    end;
end;
typecode = 'SBML_MODIFIER_SPECIES_REFERENCE';

bSBML = 0;

% check that Model is a structure
bSBML = isstruct(SBMLStructure);

% check it contains each of the fields listed
index = 1;
while (bSBML == 1 && index <= nNumberFields)
    bSBML = isfield(SBMLStructure, SBMLfieldnames(index));
    index = index + 1;
end;

% check that it contains only the fields listed
if (bSBML == 1)
    names = fieldnames(SBMLStructure);
    [m,n] = size(names);
    if (m ~= nNumberFields)
        bSBML = 0;
    end;
end;

% check that the typecode is correct
if (bSBML == 1)
    type = SBMLStructure.typecode;
    k = strcmp(type, typecode);
    if (k ~= 1)
        bSBML = 0;
    end;
end;
      
y = bSBML;
