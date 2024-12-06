nOn = input.nScansOn;
nOff = input.nScansOff;
nCells = size(data_dfof,2);

ntrials = size(input.tGratingDirectionDeg,2);
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
Dirs = unique(Dir);
nDirs = length(Dirs);

nori = length(Dirs)/2;
Ori = Dir;
for iori = 1:nori
    ind = find(Dir == Dirs(iori+nori));
    Ori(ind) = Dirs(iori);
end
Oris = unique(Ori);

%%
% Oris(1) = 360;
% 
% for iori = 1:length(Ori)
%     if Ori(iori) == 0
%         Ori(iori) = 360;
%     end
% end

%%
for n = 1:nCells
resp_win = 70:90;

% i = my_ind(1);

a = squeeze(mean(data_dfof(resp_win,:,:),1));
b = a(i,:);

% %only resp cells
% aaa = a(good_ind,:);
% bbb = aaa(i,ind);

my_x_vals = [];
my_y_vals = [];


for iori = 1:length(Oris)

    ind = find(Ori == Oris(iori));
    % indices_master(ind) = ind;
    b_ind = b(ind);   
    for trial = 1:length(ind)
        my_x_vals(iori,trial) = (b_ind(trial))*(cos(2*deg2rad(Oris(iori))));
        my_y_vals(iori,trial) = (b_ind(trial))*(sin(2*deg2rad(Oris(iori))));
    end

end

x_final = (sum(sum(my_x_vals)))/(ntrials);
y_final = (sum(sum(my_y_vals)))/(ntrials);
po_rads = (atan(y_final/x_final))/2;
po_deg(n) = rad2deg(po_rads);
%po_deg = rad2deg(po_rads);
end
%po_deg(find(po_deg<0)) = 180+po_deg(find(po_deg<0));

%%
figure; scatter(prefOri(1,:),po_deg)
xlabel('von Mises')
ylabel('Vector Sum')
%%

my_ind = find(prefOri(1,:)>55 & prefOri(1,:)<120);
weird_pos = po_deg(my_ind);

figure; scatter(prefOri(1,my_ind),weird_pos)
xlabel('von Mises')
ylabel('Vector Sum')
%%

figure;
plot(a(3,:))


%%
%nicholas solution

t = 0:22.5:167.5;
x = repmat(t, [1 10]); %standin for ori on each trial- 80 trials
y = rand(100,80); %standin for response of 100 neurons on each trial
theta = zeros(1,100);
for ii = 1:nCells
z= sum(y(ii,:).*cos(2.*deg2rad(x))+y(ii,:).*i.*sin(2.*deg2rad(x)));
theta(ii) = 0.5*rad2deg(angle(z));
end

%%

x = Ori;
y = a;
for ii = 1:nCells
z= sum(y(ii,:).*cos(2.*deg2rad(x))+y(ii,:).*i.*sin(2.*deg2rad(x)));
theta(ii) = 0.5*rad2deg(angle(z));
end

theta(find(theta<0)) = 180+theta(find(theta<0));


%%

figure;
scatter(prefOri(1,:),theta)