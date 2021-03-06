
filelist='exp_list_scmr7_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\colormaps\scmr7\';
InDirectory = 'E:\users\kenichi\scmr7\pixelmaps\'

for i=1:length(s)
    if (isempty(strfind(s(i).stim_code,'Img_flash')) & isempty(strfind(s(i).stim_code,'sweepbar')))
        continue;
    end
    load(fullfile(InDirectory,s(i).filename,'dir_dF_sm.mat'));
    colormap=getColormap(dir_dF_sm);
    imwrite(colormap, fullfile(outDirectory, [s(i).filename, '.bmp']));
end
