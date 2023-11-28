expt(1).mouse = 'i2052';
expt(1).day(1).date = '220201';
expt(1).day(2).date = '220203';
expt(1).day(1).run = '003';
expt(1).day(2).run = '004';

expt(2).mouse = 'i2053';
expt(2).day(1).date = '220208';
expt(2).day(2).date = '220210';
expt(2).day(1).run = '006';
expt(2).day(2).run = '002';

expt(3).mouse = 'i2058';
expt(3).day(1).date = '220330';
expt(3).day(2).date = '220401';
expt(3).day(1).run = '003';
expt(3).day(2).run = '002';

expt(4).mouse = 'i2062';
expt(4).day(1).date = '220427';
expt(4).day(2).date = '220429';
expt(4).day(1).run = '002';
expt(4).day(2).run = '002';

expt(5).mouse = 'i2067';
expt(5).day(1).date = '220625';
expt(5).day(2).date = '220627';
expt(5).day(1).run = '004';
expt(5).day(2).run = '002';

expt(6).mouse = 'i2070';
expt(6).day(1).date = '220802';
expt(6).day(2).date = '220804';
expt(6).day(1).run = '003';
expt(6).day(2).run = '002';

%%
fnbase = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\ACh','Data','2p_data');
fnX = ['X:\\All_staff\home\ACh\Data\2p_data'];

for iexp = 2
    mouse = expt(iexp).mouse;
    disp(mouse)
    for id = 1
        date = expt(iexp).day(id).date;
        run = expt(iexp).day(id).run;
        load(fullfile(fnX,mouse,date,run,[run '_000_000_eye.mat']))
        fnout = fullfile(fnbase,mouse,date,run);
        mkdir(fnout)
        save(fullfile(fnout,[run '_000_000_eye.mat']),'data','abstime','time')        
    end
end






