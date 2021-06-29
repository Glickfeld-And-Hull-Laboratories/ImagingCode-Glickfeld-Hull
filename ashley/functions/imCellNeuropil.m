function f = imCellNeuropil(cellmask,buffersize,neuropilsize)
% creates a neuropil mask (FOV) around each selected cell 
%       cellmask == FOV of cell ROIs
%       buffersize = width of ring around each cell
%       neuropilsize = width of neuropil ring

seBR = strel('disk', buffersize); %buffer ring
seNPR = strel('disk', buffersize+neuropilsize); %neuropil ring

siz = size(cellmask);
neuropil = zeros(siz(1),siz(2),length(unique(cellmask))-1);
for i = 1:(length(unique(cellmask))-1)
    tempmask = cellmask == i; %logical
    tempbufferring = cellmask == i;
    tempNPcells = cellmask;
    tempNPcells(tempNPcells == i) = 0; %dilate other cells so that you can subtract overlap
    cellBuffer = imdilate(tempNPcells,seBR);
    tempmask = imdilate(tempmask,seNPR) - imdilate(tempbufferring,seBR);
    tempmask = tempmask & ~cellBuffer;
    tempmask(find(isnan(tempmask))) = 0;
    neuropil(:,:,i) = tempmask;
end

cellsMaskCell = unique(cellmask);
NP = ones(size(neuropil));
for i = 1:(length(unique(cellmask))-1)
    NP(:,:,i) = neuropil(:,:,i)*i;
end
cellsNeuropil = unique(NP);
cellsNeuropil(1) = [];

missingcells = setdiff(cellsMaskCell,[0; cellsNeuropil]);

if ~isempty(missingcells)
    tempmask = zeros(size(cellmask));
    for i = missingcells
        neuropil(:,:,i) = max(neuropil(:,:,cellsNeuropil),[],3);
    end
    disp('Cells are missing neuropil, added all other neuropil')
end

f = neuropil;

end


