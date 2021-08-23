fid = fopen('RTsAdj.txt');
RTsAdj = fscanf(fid, '%f');
fclose(fid);

fid = fopen('RotEncod7Adj.txt');
RotEncod7Adj = fscanf(fid, '%f');
fclose(fid);

fid = fopen('NoJuiceAdj.txt');
NoJuiceAdj = fscanf(fid, '%f');
fclose(fid);

fid = fopen('LaserStimAdj.txt');
LaserStimAdj = fscanf(fid, '%f');
fclose(fid);
 
fid = fopen('FirstLicksEpochsAdj.txt');
FirstLicksEpochsAdj = fscanf(fid, '%f');
fclose(fid);

fid = fopen('FirstJuiceAdj.txt');
FirstJuiceAdj = fscanf(fid, '%f');
fclose(fid);

fid = fopen('AllLicksAdj.txt');
AllLicksAdj = fscanf(fid, '%f');
fclose(fid);