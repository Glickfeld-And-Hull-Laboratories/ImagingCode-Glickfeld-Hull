function out = findIsilon;

if ispc
    out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\';
elseif isunix
    out = '~/GlickfeldLabShare/All_Staff/';
else
    error('not windows or linux')
end
    