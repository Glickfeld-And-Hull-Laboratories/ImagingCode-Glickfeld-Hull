% This file contains some basic startup code for analyzing The CC dataset.
% It navigates to the folder, loads the data, and add paths where basic
% reused functions are stored. Then, you can run your blocks of code to do
% analysis and create figures.

% The location & filename of the saved data.
DataLocation = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\marie\CCconcat\MeanSubtract_2412';
DataName = 'EssentialWorkspace2.mat';

% Navigate to the location of the Workspace with the saved data & load it.
cd(DataLocation);
load(DataName);

% this code adds paths where I keep code. If the code is located somewhere
% else, change the paths.
addpath(genpath('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\marie\originalCode'));
addpath(genpath('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\marie\Kilosort\Kilosort2-noCAR'));
