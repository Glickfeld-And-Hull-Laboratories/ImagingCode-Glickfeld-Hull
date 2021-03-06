function [expmodel] = SBmergemodexp(model, experiment)
% SBmergemodexp: combines an experiment with a model, and produces a "merged 
% model", as output. The original model and experiment arguments can be
% given as textfiles, or SBexperiment objects. The output model is an
% SBmodel.
%
% DESCRIPTIONS + SYNTAX:
% ======================
% To fill in!
%
% USAGE:
% ======
% [expmodel] = SBmergemodexp(model, experiment)        
% [expmodel] = SBmergemodexp(model, experimentfile)        
%
% model: SBmodel 
% experiment: An SBexperiment object describing an experiment that should be done
%             with the model
% experimentfile: String with the name of an experiment file
%
% Output Arguments:
% =================
% Merged model containing the original model and the experimental settings.

% Information:
% ============
% Copyright (C) 2005-2007 Fraunhofer Chalmers Centre, Gothenburg, Sweden
% Main author: Henning Schmidt
% 
% Changes for the SBTOOLBOX2:
% 1/1/2008  Henning Schmidt, henning@sbtoolbox2.org
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


expmodel = [];
time = 0;   % per default time=0 is assumed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input arguments 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isSBmodel(model),
    error('The first input argument needs to be an SBmodel.');
else
    modelstruct = SBstruct(model);
end
if ~isSBexperiment(experiment),
    error('The first input argument needs to be an SBexperiment.');
else
    expstruct = SBstruct(experiment);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the output model structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newmodstruct = modelstruct;                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all initial conditions and parameters in the workspace. Do not
% define VARIABLES and REACTIONS!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kkkloopnotuseinmodel = 1:length(newmodstruct.states),
    if ~isempty(newmodstruct.states(kkkloopnotuseinmodel).initialCondition),
        eval([newmodstruct.states(kkkloopnotuseinmodel).name '=' num2str(newmodstruct.states(kkkloopnotuseinmodel).initialCondition) ';']);
    else
        error('Initial condition for state ''%s'' is undefined.',newmodstruct.states(kkkloopnotuseinmodel).name);
    end
end
for kkkloopnotuseinmodel = 1:length(newmodstruct.parameters),
    if ~isempty(newmodstruct.parameters(kkkloopnotuseinmodel).value),
        eval([newmodstruct.parameters(kkkloopnotuseinmodel).name '=' num2str(newmodstruct.parameters(kkkloopnotuseinmodel).value) ';']);
    else
        error('Value for parameter ''%s'' is undefined.',newmodstruct.parameters(kkkloopnotuseinmodel).name);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now just evaluate the paramicsettings sequentially 
% and add the new models to the structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kkkloopnotuseinmodel = 1:length(expstruct.paramicsettings)
    if ~isempty(expstruct.paramicsettings(kkkloopnotuseinmodel).formula),
        try
            % get value and assign value to variable in workspace
            valuenotuseinmodel = eval(expstruct.paramicsettings(kkkloopnotuseinmodel).formula);
            eval([expstruct.paramicsettings(kkkloopnotuseinmodel).name '=' num2str(valuenotuseinmodel) ';']);
        catch
            if expstruct.paramicsettings(kkkloopnotuseinmodel).icflag==1,
                text = sprintf('Error in initial condition setting for state ''%s'':\n\n"%s"',expstruct.paramicsettings(kkkloopnotuseinmodel).name,lasterr);
                text = sprintf('%s\n\nNote: Only a models states and parameters can be used in mathematical\nexpressions for the initial parameter and state settings.',text);
                error(text);
            else
                text = sprintf('Error in initial condition setting for parameter ''%s'':\n\n"%s"',expstruct.paramicsettings(kkkloopnotuseinmodel).name,lasterr);
                text = sprintf('%s\n\nOnly a models states and parameters can be used in mathematical\nexpressions for the initial parameter and state settings.',text);
                error(text);
            end
        end
    else
        error('Formula for initial condition for state ''%s'' is undefined: %s',expstruct.paramicsettings(kkkloopnotuseinmodel).name,lasterr);
    end
    % value determined (valuenotuseinmodel). Add it to the model structure
    if expstruct.paramicsettings(kkkloopnotuseinmodel).icflag==1,
        % if initial condition then search states
        indexnotuseinmodel = strmatch(expstruct.paramicsettings(kkkloopnotuseinmodel).name,{newmodstruct.states.name},'exact');
        if isempty(indexnotuseinmodel),
            error('Initial condition for ''%s'' defined in experiment but state does not exist in the model.',expstruct.paramicsettings(kkkloopnotuseinmodel).name);
        end 
        newmodstruct.states(indexnotuseinmodel).initialCondition = valuenotuseinmodel;
        newmodstruct.states(indexnotuseinmodel).notes = expstruct.paramicsettings(kkkloopnotuseinmodel).notes;
    else
        % if not initial condition then search in parameters 
        indexnotuseinmodel = strmatch(expstruct.paramicsettings(kkkloopnotuseinmodel).name,{newmodstruct.parameters.name},'exact');
        if ~isempty(indexnotuseinmodel),
            % only update value if parameter appears in the model (help
            % variables are handled fine).
            newmodstruct.parameters(indexnotuseinmodel).value = valuenotuseinmodel;
            newmodstruct.parameters(indexnotuseinmodel).notes = expstruct.paramicsettings(kkkloopnotuseinmodel).notes;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the parameter changes to the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copy old variables to be able to add the new variables (parameters) first
oldvariables = newmodstruct.variables;
newmodstruct.variables = [];
% Add new variables as first in the list
for k = 1:length(expstruct.parameterchanges),
    % check that no parameter is defined twice
    if ~isempty(strmatch(expstruct.parameterchanges(k).name,{expstruct.paramicsettings.name},'exact')),
        error('Parameter ''%s'' is defined twice (settings and changes).',expstruct.parameterchanges(k).name);
    end
    % delete this parameter from the model parameters
    index = strmatch(expstruct.parameterchanges(k).name,{newmodstruct.parameters.name},'exact');
    if isempty(index),
        error('Parameter ''%s'' not present in the model but changed in the experiment.',expstruct.parameterchanges(k).name);
    end
    newmodstruct.parameters(index) = [];
    newmodstruct.variables(end+1).name = expstruct.parameterchanges(k).name;
    newmodstruct.variables(end).formula = expstruct.parameterchanges(k).formula;
    newmodstruct.variables(end).notes = expstruct.parameterchanges(k).notes;
    newmodstruct.variables(end).type = '';
    newmodstruct.variables(end).compartment = '';
    newmodstruct.variables(end).unittype = '';
end
% Add old variables at the end
for k = 1:length(oldvariables),
    newmodstruct.variables(end+1).name = oldvariables(k).name;
    newmodstruct.variables(end).formula = oldvariables(k).formula;
    newmodstruct.variables(end).notes = oldvariables(k).notes;
    newmodstruct.variables(end).type = oldvariables(k).type;
    newmodstruct.variables(end).compartment = oldvariables(k).compartment;
    newmodstruct.variables(end).unittype = oldvariables(k).unittype;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add experiment events to the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(expstruct.stateevents)
    newmodstruct.events(end+1) = expstruct.stateevents(k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return the experiment model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expmodel = SBmodel(newmodstruct);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The "find index" function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [index] = find_index(name_array, comp_string)
k = 0;
l = length(name_array);
index = inf;
while k < l
    k = k+1;
    if strcmp(comp_string, name_array(k))
        index = k;
    end
end
if index == inf,
    error('The experiment description contains the element ''%s'' that is not present in the model.',comp_string);
end
return