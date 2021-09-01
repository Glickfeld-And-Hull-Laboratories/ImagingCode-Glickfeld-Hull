mouse = 'TH09'; %enter the mouse number
date = '210223'; %enter the date
user = ''; %enter your name here


base = 'D:2pdata\'; %will change this to be the 2p data folder

function [outputArg1,outputArg2] = untitled3(inputArg1,inputArg2);

data_path = fullfile(base,user,mouse,date,run);

cd(data_path)
load([run '.mat'])
data = sbxread([run],0, nframes);

data_r = squeeze(data(2,:,:,:));
data_r_avg = mean(data_r(:,:,:),3);

%register the green channel
[out data_r_reg] = stackRegister(data_g,data_g_avg);
data_r_reg_avg = mean(data_r_reg,3);


figure
imagesc(data_g_reg_avg);
title(strcat(mouse,'_',run,'_red (registered)');
colormap grey
colorbar
print(fullfile(out_path, [date '_' mouse '_runs-' run '_red_greyscale.pdf']),'-dpdf','-bestfit')

end

