 function[]  = quickGreenFOV(run, nframes);
    %user will enter the run number (all 9 digits) and the number of
    %frames.
    

    base = 'Z:\home\';
    %base = 'D:2pdata\'; %will change this to be the 2p data folder

    %data_path = fullfile(base,user,mouse,date,run);
    data_path = fullfile(base,user,'data',mouse,'2P', date, run);

    cd(data_path)
    load([run '_000_001.mat'])
    data = sbxread([run '_000_001'],0, nframes);

    data_g = squeeze(data(1,:,:,:)); %get the green channel
    data_g_avg = mean(data_g(:,:,:),3); %average it
    %register the green channel
    [out data_g_reg] = stackRegister(data_g,data_g_avg);
    data_g_reg_avg = mean(data_g_reg,3);
    
    figure
    imagesc(data_g_reg_avg);
    title(strcat(mouse," "',run," green"));
    colormap gray
    colorbar
    %save in current directory
    print(fullfile(data_path, [date '_' mouse '_run_' run '_green_greyscale.pdf']),'-dpdf','-bestfit')

end


