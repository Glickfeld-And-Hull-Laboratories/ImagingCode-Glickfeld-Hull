
function HistoOfflineAnalysis_SB3
%
%COMMAND LINE FOR VIEWING:
%
%figure(1); imagesc(log10(sum(cat(3,OUTPUT.Bkg(:).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('background');
%figure(2); imagesc(log10(sum(cat(3,OUTPUT.Hemisphere(1).Block(:,:,:).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('Left Hemisphere');
%figure(3); imagesc(log10(sum(cat(3,OUTPUT.Hemisphere(2).Block(:,:,:).count2D),3))); colorbar; set(gca,'ydir','normal'); ylabel('Virus'); xlabel('dye'); title('Right Hemisphere');


%set parameters
Mag = 1;  %=1 is the full spatial resolution;
virus = 2; %dTomato channel = viral expression
dye   = 1; %647 channel = dye capture

%spatially divide into NxM grid
NROW = 3;  %dorsal to ventral
NCOL = 4;  %medial to lateral


SAVEFILE = [pwd '\HistoIX83_Analysis_SB03.mat'];
% if isfile(SAVEFILE)
%      disp('Analysis already complete');
%      return;
% end

if ~isfile('HistoIX83.mat')
    HistoIX83_Merge;
end

%prior analysis done
HISTO = load('HistoIX83.mat');
HISTO = HISTO(1).HISTO;
HISTO.RootDir = pwd;

TotalSection = 0;
for k=1:length(HISTO.SubDir)
    TotalSection = TotalSection + length(HISTO.SubDir(k).RoiCell);
end
    
OUTPUT = [];

%Initialize Histogram
OUTPUT.ch           = [virus dye];
gamma=10; E = 65536*((0:500)/500).^(gamma); E=unique(round(E));  %generates ~300 bins, with high resolution for low intensity (gamma = 10)
OUTPUT.binEdges        = E;  %intensity resolution
OUTPUT.binIntensity = 0.5*OUTPUT.binEdges(1:end-1) +  0.5*OUTPUT.binEdges(2:end);  %intensity of each bin.
OUTPUT.NROW = NROW;
OUTPUT.NCOL = NCOL;


%dock new figure windows
set(0,'DefaultFigureWindowStyle','docked')

nSection = 0;
for k=1:length(HISTO.SubDir)  %loop through each brain section
    
    SectionROI = length(HISTO.SubDir(k).RoiCell);  %number of ROI pairs for this section (can be 0, or 1 pair;  unlikely to be more than 1 pair)
    if SectionROI>0
        %load the raw data at full resolotion
        %NOTE: bfopen only works if HISTO was saved on the computer you want to use (so rootdir is the same)
        [imageData] = bfopenIX83([HISTO.RootDir '\' HISTO.SubDir(k).Name], Mag);
        try
            if HISTO.bVS200
                imageData = rot90(imageData, -1); %counterclockwize by 90
            end
        catch
        end
        nCol = size(imageData,2);
               
        %show the image (DAPI channel, which is alwyas last), just for purposes of visual reassurance
        close all;
        imagesc(imageData(:,:,end));
        set(gca,'ColorMap', gray);
        set(gca,'ydir','reverse'); axis image; hold on;
        title(['Sample #'  num2str(k)  '  of  '  num2str(length(HISTO.SubDir))]); %  '  ' HISTO.SubDir(k).Name]);
        drawnow;
        
        BW_BKG = [];
    end
    
    for j=1:SectionROI  %loop through the ROI pairs in this sample
        nSection = nSection + 1;
        disp(['Analyzing Image ' num2str(k) ' of ' num2str(length(HISTO.SubDir)) '.  Section#' num2str(nSection) ' of ' num2str(TotalSection)]);
        
        Fields = {'Nuc'  'Cyto'}; %These names are a historical artifact.  Here, LEFT  = Nuc;  RIGHT = Cyto
        for nHemisphere=1:length(Fields)
            ROI = HISTO.SubDir(k).RoiCell(j).(Fields{nHemisphere});  %first LEFT, then RIGHT
            ROI.Pos = ROI.Pos/(2^Mag);  %convert to current magnification
            
            %draw the ROI and get the ROI mask (pixels inside)
            tmpH = feval(ROI.Type, gca, ROI.Pos);
            drawnow;
            BW = createMask(tmpH);
            drawnow;
            if isempty(BW_BKG)
                BW_BKG = (~BW);  %pixels not in this ROI
            else
                BW_BKG = (~BW) & BW_BKG;  %pixels not in any ROI so far
            end
            
            %section up the ROI into NROW x NCOL blocks
            TMP = find(sum(BW,1)); 
            ColEdges = round(min(TMP):((1+max(TMP)-min(TMP))/NCOL):1+max(TMP));
            TMP = find(sum(BW,2));
            RowEdges = round(min(TMP):((1+max(TMP)-min(TMP))/NROW):1+max(TMP));

            %determine if we are on the LHS or RHS of image
            %show symbol on dorso-medial corner.
            if ColEdges(1)-1 > nCol-ColEdges(end)
                bFlipCol = 0;  %we are on the right side, so naturally going from medial to lateral
                plot(ColEdges(1), RowEdges(1), 'ow', 'linewidth', 2);
            else
                bFlipCol = 1;  %we are on the left side, so need to flip columns to go from medial to lateral
                plot(ColEdges(end), RowEdges(1), 'ow', 'linewidth', 2);
            end

            %plot grid
            for nC = 1:length(ColEdges)
                plot(ColEdges([nC nC]), RowEdges([1 end]), ':w');
            end
            for nR = 1:length(RowEdges)
                plot(ColEdges([1 end]), RowEdges([nR nR]), ':w');
            end
            
            %Get the fancy new curvy grid
            drawnow;
            ROIset = ROI_Subdivide(ROI, NROW, NCOL);
            if bFlipCol
                %left of image, need to change columns order to go from medial-->lateral
                ROIset = ROIset(:, end:-1:1);
            end

            for nR = 1:length(RowEdges)-1
                for nC = 1:length(ColEdges)-1
                    
                    tmpH = feval(ROIset(nR,nC).Type, gca, ROIset(nR,nC).Pos);
                    setColor(tmpH,'r')
                    drawnow;
                    BW2 = createMask(tmpH);
                    drawnow;
                    %plot(ROIset(nR,nC).Pos(:,1), ROIset(nR,nC).Pos(:,2), '.r');
                    %drawnow;

                    %pixels in curvy grid (also make sure they are in the original ROI)
                    I = find( BW & BW2  ); 
                    
                    Pixels = zeros(length(I), 2);
                    for q = 1:length(OUTPUT.ch) %for the channels we care about
                        tmp = imageData(:,:,OUTPUT.ch(q)); %full image, single color
                        Pixels(:,q) = tmp(I);
                    end
                    
                    OUTPUT.Hemisphere(nHemisphere).Block(nR, nC, nSection).ROI = ROIset(nR,nC);
                    OUTPUT.Hemisphere(nHemisphere).Block(nR, nC, nSection).count2D = histcounts2(Pixels(:,1), Pixels(:,2), OUTPUT.binEdges, OUTPUT.binEdges);
                end
            end
        end  %finished one Section/Hemisphere
        
        %BKG histogram for this section:
        I = find(BW_BKG);
        
        Pixels = zeros(length(I), 2);
        for q = 1:length(OUTPUT.ch) %for the channels we care about
            tmp = imageData(:,:,OUTPUT.ch(q)); %full image, single color
            Pixels(:,q) = tmp(I);
        end        
        OUTPUT.Bkg(nSection).count2D = histcounts2(Pixels(:,1), Pixels(:,2), OUTPUT.binEdges, OUTPUT.binEdges);
        
    end
end
save(SAVEFILE, 'OUTPUT', '-v7.3');







