function status = save_all_open_figs(path)
% status = save_all_open_figs(path,identifier)
%   Saves all open figures to a given directory as pdf.
%   PATH should be a file path provided as characters.
%   IDENTIFIER should be numerical 1 or 2. 1 will result in figures saved
%   by figure NAME, 2 will result in figures saved by figure TITLE.

figHandles = findall(0,'Type','figure'); 
this_fig = get(figHandles,'Name');

for i = 1:numel(figHandles)
    t = this_fig{i,1};
    saveas(figHandles(i),fullfile(path,[t '.pdf']));
end

status = "Figures saved!";

end