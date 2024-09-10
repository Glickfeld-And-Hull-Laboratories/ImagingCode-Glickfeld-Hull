function experiment_sessions = query_expt(mouse)
% experiment_sessions = query_expt(mouse_ID)
%   MOUSE_ID should be the mouse number (without the "i")
%   EXPERIMENT_SESSIONS returns all sessions associated with the mouse in a
%   double array column


ds = 'DART_expt_info'; 
eval(ds);
mouse_id = convertCharsToStrings(['i' num2str(mouse)]);
subnums = vertcat(expt.SubNum);
nLoops = size(subnums,1);
substr = strings(nLoops,1);
for row = 1:nLoops
    substr(row,1) = convertCharsToStrings(subnums(row,:));
end
isThere = strcmp(substr, mouse_id);
if sum(isThere) == 0
    error('Mouse ID does not exist. Enter a valid mouse ID.')
else
    experiment_sessions = find(isThere);
end

end