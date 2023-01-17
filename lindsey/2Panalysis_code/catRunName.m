function run_str = catRunName(ImgFolder, nrun);
%creates name string for saving files
if ~iscell(ImgFolder)
    ImgFolder = mat2cell(ImgFolder,1);
end
run_str = ['runs-' ImgFolder{1}];
if nrun>1
    for irun = 2:nrun
        run_str = [run_str '-' ImgFolder{irun}];
    end
end