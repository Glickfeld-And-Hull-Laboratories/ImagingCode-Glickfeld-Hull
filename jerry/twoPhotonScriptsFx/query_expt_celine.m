function experiment_sessions = query_expt_celine(mouse)
% experiment_sessions = query_expt(mouse_ID)
%   For Jerry's DART experiment metadata. Returns corresponding DART experiment session numbers for
%   the current mouse. 
%   MOUSE_ID should be the mouse number (without the "i")
%   EXPERIMENT_SESSIONS returns all sessions associated with the mouse in a
%   double array column


ds = 'DART_V1_contrast_ori_Celine'; 
eval(ds);
mouse_id = convertCharsToStrings(['i' num2str(mouse)]);
% subnums = vertcat(expt.SubNum);
subnums = cell(length(expt),1);
for i = 1:length(expt)
    subnums{i} = expt(i).SubNum;
end

nLoops = size(subnums,1);
substr = strings(nLoops,1);
for row = 1:nLoops
    if isempty(subnums{row})
        substr(row,1) = 'empty';
    else
        substr(row,1) = convertCharsToStrings(subnums{row});
    end
end

isThere = strcmp(substr, mouse_id);

if sum(isThere) == 0
    error(['Mouse ID does not exist. Enter a valid mouse ID. Mouse ID should be the 4-digit number without the "i".'])
else
    experiment_sessions = find(isThere);
end

end