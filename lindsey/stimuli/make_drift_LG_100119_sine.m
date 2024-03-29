%make a square wave grating, save it.. 


LUM_HI = .8; %.95; %set to 0 and 1 for prev version
LUM_LO = 0; %.4;
LUM_MEAN = .3; %figure this out impirically using calibration
                  %using full field blanks (see end of this file,
                  %after the word 'break', 
                  %for how to make these)
		  
scale1 = 4; %scales the pixelation so that movies are smaller and
            %run better. this is good to do for square wave
            %movies.. 
	    
%convert degrees into pixels: 
%at 20cm, 1cycle = 20 degrees => tan(20)=O/20
%so O = 20*tan(20*pi/180) %in cm
pix_per_cm = 24./scale1;

Ly = 253; %choose appropriate values for your monitor
Lx = 321;

PWD_visfiles = ['G:\users\lindsey\svn\users\lindsey\stimuli'];

a = dir(PWD_visfiles);
if isempty(a)
  eval(['!mkdir ',PWD_visfiles]);
end

ori_vec = [0]*pi/180; %orientation of gratings
Dist = 20; %cm from screen 
f_vec = [.05]; %cycles per degree [.025 .05 .1]
Contrast = 100; %percent

color = 'w'; %can be 'b' or 'w' right now.. 
	     %now make a 3color version: 
if strcmp(color,'w')
  Cmap1 = [1 1 1];
elseif strcmp(color,'b');
  Cmap1 = [0 0 1];
end
 
DriftRateMat = 1; % spatial frequency

phase = 0*pi/180; %phase of grating (in rad)

X = [1:Lx]'*ones(1,Ly);
Y = ones(Lx,1)*[1:Ly];
for count_f = 1:length(f_vec)
  f = f_vec(count_f);
  f_txt = num2str(1/f);
  
%  Opp =  Dist*tan((1/f)*pi/180)* pix_per_cm %length of 1 cycle, in pixels
  Opp =  Dist*tan(1/(2*f)*pi/180)* pix_per_cm %length of 1 cycle, in pixels

% generate gabor filter
sigma = 42;
gamma = .79;
theta = ori_vec;
psi = 1;
lambda = 1;

sigma_x = sigma;
sigma_y = sigma/gamma;
 
% Bounding box
nstds = 3;
xmax = max(abs(nstds*sigma_x*cos(theta)),abs(nstds*sigma_y*sin(theta)));
xmax = ceil(max(1,xmax));
ymax = max(abs(nstds*sigma_x*sin(theta)),abs(nstds*sigma_y*cos(theta)));
ymax = ceil(max(1,ymax));
xmin = -xmax; ymin = -ymax;
[x,y] = meshgrid(xmin:xmax,ymin:ymax);
 
% Rotation 
x_theta=x*cos(theta)+y*sin(theta);
y_theta=-x*sin(theta)+y*cos(theta);


gb= 1/(2*pi*sigma_x *sigma_y) * exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);
gb_gradient = LUM_HI.*(gb./max(max(gb)));
  %first, generate a sine wave grating.. 
  Nphase = ceil(2./f);
  Big_Z3 = zeros(Nphase,Lx,Ly,3);
  for count_phase = 1:Nphase
      phase = count_phase./Nphase*2*pi;
      for count_ori = 1:length(ori_vec);
        ori = ori_vec(count_ori);
        ori_txt = num2str(round(ori*180/pi*100));
        if length(ori_txt) == 1
            ori_txt = ['00',ori_txt];
        elseif length(ori_txt) == 2
            ori_txt = ['0',ori_txt];
        elseif length(ori_txt) == 3
            ori_txt = [ori_txt];
        end
        Z = sin((X*sin(ori) + Y*cos(ori))./(2*Opp)*2*pi + phase);
    %  imagesc(Z2); colormap('gray'); axis image; pause(.2)
    %Z2_use = (Z2a-.5)*Contrast/100 + .5;
    %assume Z2a is all 0's and 1's
    Z2a = (Z+1)./2;
    Z2_use = zeros(size(Z));
    Z2_use = Z2a.*(gb_gradient);
    Z3 = zeros(Lx,Ly,3);
    for count = 1:3
      Z3(:,:,count) = Z2_use*Cmap1(count);
    end
    image(Z3); axis image; pause(.2)
    Big_Z3(count_phase,:,:,:) = Z3;
      end
  end
  
  
  %now make an avi
  Ncyc = f*80;
  Nframes = Ncyc * Nphase;
  avi_Z3 = zeros(Nframes,Lx,Ly,3);
  
  for count_drift = 1:length(DriftRateMat)
    DriftRate = DriftRateMat(count_drift);
    count_phase2a = 0;
    for n0 = 1:Ncyc
      for count_phase = 1:Nphase
	count_phase2a = count_phase2a + 1;
	
	if DriftRate == 1
	  count_phase2 = rem(count_phase2a,Nphase) + 1;
	elseif DriftRate == 0
	  count_phase2 = 1;
	elseif DriftRate == 2
	  %slow down by a factor of 2: 
	  count_phase2 = rem(ceil(count_phase2a/2),Nphase) + 1
	elseif DriftRate == 3
	  count_phase2 = rem(ceil(count_phase2a/3),Nphase) + 1
%	  count_phase2 = count_phase;
	end
	
	
	Z3B = squeeze(Big_Z3(count_phase2,:,:,:));
	h = figure(gcf); clf;
	
	h2 = image(Z3B); axis image;
	axis off;
	pause(.1);
	%M((n0-1)*2 + count_rev) = getframe(gcf,[440 314 560 420]);;
	M((n0-1)*Nphase + count_phase) = im2frame(Z3B);
      end 
    end
 
  file_save0 = [PWD_visfiles,'100119_vid_drift_sinewavegrat_', ori_txt,'ori_',color,'_',...
      f_txt,'SF_',num2str(DriftRate),'TF_',num2str(Contrast),'cont_',...
      num2str(LUM_LO),'_',num2str(LUM_HI),'LUM','.avi'];
  speed = Nphase*DriftRateMat;
  movie2avi(M,file_save0,'compression','none');
  end
end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


break

%make a blank screen with the medium color; 

mean_lum = .7565; %.7565 from June20.7; %.35;
Z4 = zeros(Lx,Ly,3);
for count = 1:3
  Z4(:,:,count) = mean_lum*Cmap1(count);
end

strP = ['_Ppt',num2str(mean_lum*100),'_'];


file_save = [PWD_visfiles,'Jmeanlum',strP,'_c',color,'.jpg'];
imwrite(Z4,file_save,'bitdepth',8,'Quality',100);

file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'.tif'];
imwrite(Z4,file_save,'compression','none');

%make an avi of teh blank screen
file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'.avi'];
M =  im2frame(Z4);
movie2avi(M,file_save,'compression','None');



%make files of many luminances: 

for mean_lum = .3:.01:1
  Z4 = zeros(Lx,Ly,3);
  for count = 1:3
    Z4(:,:,count) = mean_lum*Cmap1(count);
  end
  
  strP = ['_Ppt',num2str(mean_lum*100),'_'];
  
  file_save = [PWD_visfiles,'Jmeanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.jpg'];
  imwrite(Z4,file_save,'bitdepth',8,'Quality',100);

  file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.tif'];
imwrite(Z4,file_save,'compression','none');
  
  
  %make an avi of teh blank screen
  file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.avi'];
  M =  im2frame(Z4);
  movie2avi(M,file_save,'compression','None');
  
end


for mean_lum = .51:.001:.57; %.621:.001:.629
  Z4 = zeros(Lx,Ly,3);
  for count = 1:3
    Z4(:,:,count) = mean_lum*Cmap1(count);
  end
  
  strP = ['_Ppt',num2str(mean_lum*1000),'_'];
  
  file_save = [PWD_visfiles,'Jmeanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.jpg'];
  imwrite(Z4,file_save,'bitdepth',8,'Quality',100);

  file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.tif'];
imwrite(Z4,file_save,'compression','none');
  
  
  %make an avi of teh blank screen
  file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.avi'];
  M =  im2frame(Z4);
  movie2avi(M,file_save,'compression','None');
  
end




