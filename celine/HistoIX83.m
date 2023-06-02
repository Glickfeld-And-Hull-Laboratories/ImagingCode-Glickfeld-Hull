
function varargout = HistoIX83(varargin)
%function varargout = HistoIX83(varargin)
%
%handles time series

if nargin == 0
    LaunchGUI;
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch ERR
        rethrow(ERR);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LaunchGUI
set(0,'DefaultFigureWindowStyle','docked');
if ~isempty(whos('global','HISTO'))
    Title = 'Warning: is HistoIX83 already running?';
    Question = 'Only one instance of HistoIX83 can run at a time.  If you are sure one is not already running press OK';
    ButtonName=questdlg(Question,Title,'OK','Cancel','Cancel');
    if ~strcmpi(ButtonName, 'OK')
        return;
    end
end
%Get directory information
D = GetDirRecursive;
if isempty(D)
    errMsg = 'There are no valid .tif files in this directory, or any of its subfolders. Program will not launch.';
    waitfor(errordlg(errMsg));
    return;
end

global HISTO
HISTO = [];

%Added 4-23-15 to speed up performance.
HISTO.Mag = 4;  %smaller number = higher resolution images, but slower.  Coordinate transform is preserved
%HISTO.Mag = 3;  %

HISTO.RootDir = pwd;
HISTO.bFirstClick= true;
HISTO.bVS200 = false;

bNormalFig = 1;
if bNormalFig
    fig = figure('Name', 'HistoIX83');
    %fig = figure('Name', 'HistoIX83', 'position', [50    50   1320   1600]);
    %set(gcf, 'Pointer', 'watch'); drawnow;
else
    %create the main figure -- default position
    %fig = figure('Name', 'HistoIX83', 'NumberTitle', 'off', 'menubar', 'none', 'position', [50    50   1320   753]);
    %fig = figure('Name', 'HistoIX83', 'NumberTitle', 'off', 'position', [50    50   1320   753]);
    fig = figure('Name', 'HistoIX83', 'NumberTitle', 'off');
 HistoIX83   %set(gcf, 'Pointer', 'watch'); drawnow;
    
    %make handle visibility full so that we can create dynamic uicontrols in this function
    %then set back to 'callback' at end of initialization
    if ~bNormalFig
        set(fig, 'HandleVisibility', 'on', 'IntegerHandle', 'on');
    end
end
set(fig,'toolbar','figure');

%Set the CloseRequestFcn -- to save options on close
set(fig, ...
    'CloseRequestFcn', CallbackString('HistoIX83', '''MainFigure_CloseRequestFcn'''), ...
    'ResizeFcn', CallbackString('HistoIX83', '''MainFigure_ResizeFcn'''), ...
    'WindowButtonDownFcn', CallbackString('HistoIX83', '''MouseButtonDown'''), ...
    'WindowButtonUpFcn', CallbackString('HistoIX83', '''MouseButtonUp'''), ...
    'WindowScrollWheelFcn', @(obj, evd)  HistoIX83('MainFigure_WindowScrollWheelFcn', obj, evd), ...;
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

%Now create all the main UI objects
White = [.95 .95 .95];

HISTO.UI = []; %place holder


%now integrate this into our prior framework
N = 0;
for n=1:length(D)
    N = N+1;
    HISTO.SubDir(N).Name   = D(n).name;
    HISTO.SubDir(N).NumCh  = NaN;  %number of color channels
    HISTO.SubDir(N).Axis   =  [];  %[Xmin Xmax Ymin Ymax]
    HISTO.SubDir(N).ChData = [];  %this is a vector of data at different resolutions.  Each contains:
    %       bPlotted       binary, 1 if already rendered
    %       mag   scalar: (0 to 9)
    %       image 16bit:  (row x col x channels)
    %       absColX, absRowY: vectors indicating absolute coordinates (at highest mag)
    HISTO.SubDir(N).RoiBkg = [];
    HISTO.SubDir(N).RoiCell = [];
    HISTO.SubDir(N).bFlipLR = false;
end


%try to load data
try    
    DAT = read([HISTO.RootDir '\HistoIX83.mat']);
    
    %error-check for directory consistency.
    bErrorFileMatch = 0;
    if length(DAT.HISTO.SubDir) ~= length(HISTO.SubDir)
        bErrorFileMatch = 1;
    end
    for n=1:length(DAT.HISTO.SubDir)
        if ~strcmpi(DAT.HISTO.SubDir(n).Name, HISTO.SubDir(n).Name)
            bErrorFileMatch = 1;
        end
    end        
    if bErrorFileMatch
        %make a copy of the current file
        NewName = 'HistoIX83_copy.mat';
        while ~isempty(dir([HISTO.RootDir '\' NewName]))
            %ensure that we don't overwrite a pre-existing file
            NewName = [NewName(1:end-4) '_copy.mat'];
        end
        movefile([HISTO.RootDir '\HistoIX83.mat'],[HISTO.RootDir '\' NewName]);
        errMsg = ['Saved HistoIX83.mat has a mismatch in the number of .tif files or subdirectory structure.  Saved settings cannot be loaded; original file has been renamed to ' NewName];
        waitfor(errordlg(errMsg));
        error(errMsg);
    end
    %this code is only executed if no errors occur.
    DAT.HISTO.Mag = HISTO.Mag;
    HISTO = CopyConservedFieldValues(HISTO, DAT.HISTO, 1, []);
    HISTO.bLoadingData = 0;
    HISTO.bPlotting  = 0;
    HISTO.RootDir = pwd;  %ovverride what was saved.
catch
    
    %Standard
%     DefaultRGB = { ...
%         1   1    1   10000      0.8;
%         1   2    1   15000      0.8;
%         1   3   200   5000      0.4 };
    
    %Sasha
    DefaultRGB = { ...
        1   3    1    300      0.8;
        1   1    1   3000      0.8;
        1   2    1    500      1 };
    
    %OpenFieldBehaviorStandard
    DefaultRGB = { ...
        1   1    1   6000      0.8;
        1   2    1   15000     0.5;
        1   3    1   6000      0.8 };

    %OpenFieldBehaviorStandard with SCY5-DART
    DefaultRGB = { ...
        1   2    1     10000     1;
        1   1    10     1000     1;
        1   4    100    7000     1 };
    
    %OpenFieldBehaviorStandard with SCY5-DART
    DefaultRGB = { ...
        1   2    1     15000     0.6;
        1   1    1     10000     0.7;
        1   1    1     10000     0.7 };    
    
    %OpenFieldBehaviorStandard with AF488-DART
    DefaultRGB = { ...
        1   2    1     20000     0.7;
        1   3    1     20000     0.7;
        1   4    1     10000     1 };
    
    Fld = {'On' 'Ch' 'Min' 'Max' 'Gam'};
    for c=1:3
        for k=1:length(Fld)
            HISTO.DispRGB(c).(Fld{k}) = DefaultRGB{c,k};
        end
    end
    HISTO.bLoadingData = 0;
    HISTO.bPlotting = 0;
    HISTO.nSel = 1;
    HISTO.tSel = 1;
    HISTO.zSel = 1;  %note, this is typically overridden by dXYZ
end
%backwards compatibility -- make sure these are initialized
try
    HISTO.tSel;
catch
    HISTO.tSel=1;
end
try
    HISTO.zSel;
catch
    HISTO.zSel=1;
end

HISTO.UI.RoiMovement = 0;
HISTO.UI.Ax = axes('Units', 'pixels');
for a=1:3
    HISTO.UI.Axs(a) = axes('Units', 'pixels');
end

%HISTO.UI.SliderT = uicontrol('Style','slider','Value', 1, 'Min', 1, 'Max', 2, 'Units', 'pixels', ...
%    'Callback', CallbackString('HistoIX83','''SliderT_Callback'''), ...
%    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

%HISTO.UI.SliderZ = uicontrol('Style','slider','Value', 1, 'Min', 1, 'Max', 2, 'Units', 'pixels', ...
%    'Callback', CallbackString('HistoIX83','''SliderZ_Callback'''), ...
%    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));


HISTO.UI.ButtonSaveBackup = uicontrol('Style', 'pushbutton', 'String', 'SaveBkp', 'Units', 'pixels', ...
    'Callback', CallbackString('HistoIX83','''SaveBackup'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

Flds = {'R' 'G' 'B'};
for c=1:length(Flds)
    HISTO.UI.RGBCheck.(Flds{c}) = uicontrol('Style', 'checkbox', 'String', Flds{c}, 'Value', HISTO.DispRGB(c).On, 'Units', 'pixels', 'HorizontalAlignment', 'right', ...
        'Callback', CallbackString('HistoIX83','''RGBCheck_Callback'''), ...
        'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));
end
Flds = {'Ch' 'Min' 'Max' 'Gam'};
for k=1:length(Flds)
    HISTO.UI.Txt.(Flds{k}) = uicontrol('Style', 'text', 'String', Flds{k}, 'Units', 'pixels', 'HorizontalAlignment', 'center', ...
        'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));
    for c=1:3
        HISTO.UI.RGB(c).(Flds{k}) = uicontrol('Style', 'edit', 'String', HISTO.DispRGB(c).(Flds{k}), 'Units', 'pixels', 'HorizontalAlignment', 'center', ...
            'Callback', CallbackString('HistoIX83','''Edit_Callback'''), 'BackgroundColor', White, ...
            'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));
    end
end

HISTO.UI.ButtonREDRAW = uicontrol('Style', 'pushbutton', 'String', 'REDRAW', 'Units', 'pixels', ...
    'Callback', CallbackString('HistoIX83','''ButtonREDRAW'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

HISTO.UI.ButtonNewCellEE = uicontrol('Style', 'pushbutton', 'String', 'EE', 'Units', 'pixels', ...
    'Callback', CallbackString('HistoIX83','''ButtonNewCell_Callback''', '''imellipse''', '''imellipse'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));
HISTO.UI.ButtonNewCellEP = uicontrol('Style', 'pushbutton', 'String', 'EP', 'Units', 'pixels', ...
    'Callback', CallbackString('HistoIX83','''ButtonNewCell_Callback''', '''imellipse''', '''impoly'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));
HISTO.UI.ButtonNewCellPP = uicontrol('Style', 'pushbutton', 'String', 'PP', 'Units', 'pixels', ...
    'Callback', CallbackString('HistoIX83','''ButtonNewCell_Callback''', '''impoly''', '''impoly'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

HISTO.UI.ListBox = uicontrol('Style', 'listbox', 'String', '', 'Units', 'pixels', 'BackgroundColor', White, ...
    'Callback', CallbackString('HistoIX83','''ListBox_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

HISTO.UI.Mag = uicontrol('Style', 'popupmenu', 'String', {'Res-1' 'Res-2' 'Res-3' 'Res-4' 'Res-5'}, 'Units', 'pixels', 'BackgroundColor', White, ...
    'Value', HISTO.Mag, ...
    'Callback', CallbackString('HistoIX83','''Mag_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

HISTO.UI.bFlipLR = uicontrol('Style', 'checkbox', 'String', 'FlipLR', 'Value', 0, 'Units', 'pixels', 'HorizontalAlignment', 'right', ...
    'Callback', CallbackString('HistoIX83','''FlipLR_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

HISTO.UI.bVS200 = uicontrol('Style', 'checkbox', 'String', 'VS200', 'Value', 0, 'Units', 'pixels', 'HorizontalAlignment', 'right', ...
    'Callback', CallbackString('HistoIX83','''VS200_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  HistoIX83('MainFigure_KeyPressFcn', obj, evd));

%Resize -- setup PLUGIN_UI_RECT before loading any plugins
MainFigure_ResizeFcn;
ListBoxRefresh;
if isempty(HISTO.nSel)
    HISTO.nSel = 1;  %bug fix Feb 18, 2020
end
set(HISTO.UI.ListBox, 'Value',HISTO.nSel);

%set(gcf, 'Pointer', 'arrow'); drawnow;

%LoadData(1);

%make handle visibility full so that we can create dynamic uicontrols in this function
%then set back to 'callback' at end of initialization
%set(gcf, 'HandleVisibility', 'callback', 'IntegerHandle', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = GetDirRecursive(Path)
if nargin<1
    Path = ''; %root
end

D = []; %combined directory information -- recursive.

%first look for any *.tif files
TMP = dir([Path '*.tif']);
if ~isempty(TMP)
    %Added July 2021 -- only need to sort if some filenames are longer than
    %others.  This corrects an IX83 crash/hexadecimal naming issue.
    if length(unique(cellfun('length',{TMP.name})))>1
        %sort this folder according to the first number in each name
        for j=1:length(TMP)
            TMP(j).Num = NumbersInString(TMP(j).name);
            TMP(j).Num = TMP(j).Num(1);
        end
        [~,I] = sort([TMP.Num],'ascend');
        TMP = TMP(I);
    else
        %no sorting needed
        for j=1:length(TMP)
            TMP(j).Num = NaN;
        end
    end
    %append root directory with file name
    for j=1:length(TMP)
        TMP(j).name = [Path  TMP(j).name];
    end
    D = [D; TMP];
end

%next recurse all the subdirectories
if isempty(Path)
    DRoot = dir;
else
    DRoot = dir(Path);
end
for q=1:length(DRoot)
    if DRoot(q).isdir && isempty(findstr(DRoot(q).name, '.'))
        TMP = GetDirRecursive([Path DRoot(q).name '\']);
        D = [D; TMP];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MainFigure_ResizeFcn
%Called whenever the figure is resized

global HISTO

%Layout constants
PosMain = get(gcf,'Position');
FigW = PosMain(3);
FigH = PosMain(4);
%PosMain

S    = 5;  %space betw ax and ctrls
RGBw = 40;
ButH = 20;

BoxH = FigH - ButH - 5*ButH;
AxsH = BoxH*0.00022; %small axes height
BoxH = BoxH - AxsH*length(HISTO.UI.Axs);

BoxW = 350;  %width of list boxes in pixels
AxW = FigW-BoxW-S; 
AxH = FigH-S; %*1040/1392;

set(HISTO.UI.Ax, 'Position', [S  S  AxW  AxH], 'visible', 'off');
xlim = get(HISTO.UI.Ax, 'xlim');
ylim = get(HISTO.UI.Ax, 'ylim');
if diff(xlim)/diff(ylim) < AxW/AxH
    %xlim is too narrow -- widen it
    dX = (diff(ylim)*AxW/AxH - diff(xlim))/2;
    dX = min(dX, xlim(1));  %don't allow negative
    set(HISTO.UI.Ax, 'xlim', [xlim(1)-dX  xlim(2)+dX]);
elseif diff(ylim)/diff(xlim) < AxH/AxW
    %ylim is too narrow -- widen it
    dY = (diff(xlim)*AxH/AxW - diff(ylim))/2;
    dY = min(dY, ylim(1));  %don't allow negative
    set(HISTO.UI.Ax, 'ylim', [ylim(1)-dY  ylim(2)+dY]);
end
%axis(HISTO.UI.Ax, 'equal');


Flds = {'ButtonREDRAW'  'ButtonNewCellEE' 'ButtonNewCellEP' 'ButtonNewCellPP'};
FldW = BoxW/length(Flds);
for k=1:length(Flds)
    set(HISTO.UI.(Flds{k}), 'Position', [AxW+S+FldW*(k-1),  FigH-5*ButH, FldW, ButH]);
end
MagW = 70;
set(HISTO.UI.bVS200,           'Position', [AxW+S,        FigH-6*ButH, MagW-S, ButH]);
set(HISTO.UI.bFlipLR,          'Position', [AxW+  MagW+S, FigH-6*ButH, MagW-S, ButH]);
set(HISTO.UI.Mag,              'Position', [AxW+2*MagW+S, FigH-6*ButH, MagW-S, ButH]);
set(HISTO.UI.ButtonSaveBackup, 'Position', [AxW+3*MagW+S, FigH-6*ButH, BoxW-3*MagW, ButH]);

Flds = {'R' 'G' 'B'};
for k=1:length(Flds)
    set(HISTO.UI.RGBCheck.(Flds{k}), 'Position', [AxW+S,  FigH-ButH*(k+1), RGBw, ButH]);
end
Flds = {'Ch' 'Min' 'Max' 'Gam'};
FldW = (BoxW-RGBw)/length(Flds);
for k=1:length(Flds)
    set(HISTO.UI.Txt.(Flds{k}), 'Position', [AxW+S+RGBw+FldW*(k-1),  FigH-ButH, FldW, ButH]);
    for c=1:3
        set(HISTO.UI.RGB(c).(Flds{k}), 'Position', [AxW+S+RGBw+FldW*(k-1),  FigH-(c+1)*ButH, FldW, ButH]);
    end
end
set(HISTO.UI.ListBox, 'Position', [AxW+S, AxsH*length(HISTO.UI.Axs), BoxW, BoxH]);
for k=1:length(HISTO.UI.Axs)
    set(HISTO.UI.Axs(k), 'Position', [AxW+S, 1+(k-1)*AxsH, BoxW, AxsH], 'visible', 'off');
end
uicontrol(HISTO.UI.ListBox); %set focus to ListBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveBackup
global HISTO;
TMP = HISTO;
clear HISTO; %unlink from global variable
HISTO = TMP; %local coapy (NOT GLOBAL)
clear TMP;  %save memory
set(gcf, 'Pointer', 'watch'); drawnow;
try
    %save data without ChData
    for n=1:length(HISTO.SubDir)
        HISTO.SubDir(n).ChData = [];
        
        %analagous to DestructROI, but with no effect on actual GUI objects
        try
            HISTO.SubDir(n).RoiBkg.Pos = HISTO.SubDir(n).RoiBkg.h.getPosition;
        catch
        end
        HISTO.SubDir(n).RoiBkg = UnlinkROI(HISTO.SubDir(n).RoiBkg);
        
        for j=1:length(HISTO.SubDir(n).RoiCell)
            try
                HISTO.SubDir(n).RoiCell(j).Nuc.Pos = HISTO.SubDir(n).RoiCell(j).Nuc.h.getPosition;
            catch
            end
            HISTO.SubDir(n).RoiCell(j).Nuc = UnlinkROI(HISTO.SubDir(n).RoiCell(j).Nuc);
            
            try
                HISTO.SubDir(n).RoiCell(j).Cyto.Pos = HISTO.SubDir(n).RoiCell(j).Cyto.h.getPosition;
            catch
            end
            HISTO.SubDir(n).RoiCell(j).Cyto = UnlinkROI(HISTO.SubDir(n).RoiCell(j).Cyto);
        end
    end
    HISTO.UI = [];
    disp('saving backup ...');
    save([HISTO.RootDir '\HistoIX83_Backup.mat'], 'HISTO', '-v7.3');
    disp('done!');
catch
    Title = 'Error Saving Backup';
    errordlg('Error Saving Backup; Data not saved.',Title);
end
set(gcf, 'Pointer', 'arrow'); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MainFigure_CloseRequestFcn
global HISTO;

try
    %set(gcf, 'Pointer', 'watch'); drawnow;
    
    %save data without ChData
    for n=1:length(HISTO.SubDir)
        HISTO.SubDir(n).ChData = [];
    end
    %also clear all handles.
    ClearROIHandles;
    HISTO.UI = [];
    
    save([HISTO.RootDir '\HistoIX83.mat'], 'HISTO', '-v7.3');
    
    clear global HISTO;
    %set(gcf, 'Pointer', 'arrow'); drawnow;
    
    closereq;  %standard matlab close function
    
    %LoadDataToMatrix;
catch
    Title = 'Error closing';
    Question = 'Error closing; Data not saved.  If you are sure you want to exit press OK';
    ButtonName=questdlg(Question,Title,'OK','Cancel','Cancel');
    if ~strcmpi(ButtonName, 'OK')
        return;
    end
    %clear global HISTO;
    closereq;  %standard matlab close function
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function MainFigure_WindowScrollWheelFcn(obj, evd)
global HISTO;
return; %disable zooming for now

% axes(HISTO.UI.Ax);
% if evd.VerticalScrollCount > 0
%     zoom(1/2);
% elseif evd.VerticalScrollCount < 0
%     zoom(2);
% end
% MainFigure_ResizeFcn;
% HISTO.SubDir(HISTO.nSel).Axis = round(axis(HISTO.UI.Ax));
% ReplotAx(0,0,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MainFigure_KeyPressFcn(obj, evd)
%processes shortcut keys

global HISTO;

switch lower(evd.Key)
    case 'a'
        ButtonREDRAW;
    case 's'
        ButtonNewCell_Callback('imellipse', 'imellipse');
    case 'd'
        ButtonNewCell_Callback('impoly', 'impoly');
    case {'r' 'g'  'b'}
        SingleRGB(upper(evd.Key));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportImage

global HISTO;

[filename, pathname] = uiputfile('*.jpg', 'Export Image to JPG');
if isequal(filename,0) || isequal(pathname,0)
    return; %cancel
end

%get image
error('ExportImage code needs to be fixed');
rgb = 0; %PLACEHOLDER HISTO.SubDir(HISTO.nSel).RGB(nz).RGB;
for c=1:3
    if HISTO.DispRGB(c).On==0
        rgb(:,:,c) = 0; %black-out one color
    end
end

imwrite(rgb, [pathname '\' filename], 'jpeg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bLoadDataCalled = ReplotAx(bRedrawROI, bLoadDataIfNeeded, bDrawnow)
global HISTO;

n = HISTO.nSel;

%set FlipLR flag
set(HISTO.UI.bFlipLR, 'value', HISTO.SubDir(n).bFlipLR);
set(HISTO.UI.bVS200, 'value', HISTO.bVS200);
    
%axis(HISTO.UI.Ax, 'equal');
%axis(HISTO.UI.Ax, 'off');

if bRedrawROI
    ClearROIHandles; %clear ROI UI objects
end

%render all non-plotted images
if ~isempty(HISTO.SubDir(n).ChData)    
    %sort from low mag to high mag
    [~,I] = sort([HISTO.SubDir(n).ChData.mag]);
    HISTO.SubDir(n).ChData = HISTO.SubDir(n).ChData(I(end:-1:1));
    for k=1:length(HISTO.SubDir(n).ChData)
        if ~HISTO.SubDir(n).ChData(k).bPlotted
            
            %draw image
            axes(HISTO.UI.Ax);
            hold on;
            %on 4-23-15: scale x/y by 2^HISTO.Mag to speed up performance
            RGB = MakeRGB(HISTO.SubDir(n).ChData(k).imageData, HISTO.DispRGB);
            image(HISTO.SubDir(n).ChData(k).absColX,  HISTO.SubDir(n).ChData(k).absRowY,  RGB);
            
            %place it after annotations
            %CHILD = get(HISTO.UI.Ax, 'children');
            %CHILD = [CHILD(2:NumAnnotations+1);  CHILD(1); CHILD(NumAnnotations+2:end)];
            %set(HISTO.UI.Ax,'children', CHILD);
            
            %set flag
            HISTO.SubDir(n).ChData(k).bPlotted = 1;
            
            
            set(HISTO.UI.Ax,'ydir','reverse');
            if HISTO.SubDir(n).bFlipLR
                %'reverse'
                set(HISTO.UI.Ax,'xdir','reverse');
            else
                set(HISTO.UI.Ax,'xdir','normal');
            end
        end
    end
end

if ~bLoadDataIfNeeded 
    bLoadDataCalled = 0;
elseif ~isempty(HISTO.SubDir(n).ChData) && EnoughDataLoaded(n)
    bLoadDataCalled = 0;
else
    bLoadDataCalled = 1;
    LoadData(1);  %note, this will repetitively call ReplotAx(0,0,1)
end
    
%set(gcf, 'Pointer', 'arrow'); 
axis image
axis on;
%axis([0 150 0 150])

if bRedrawROI
    ReplotROI(bDrawnow);   %replot ROI objects
end

%determine number of ROI objects
%CHILD = get(HISTO.UI.Ax, 'children');
%NumAnnotations = 0;
%while NumAnnotations<length(CHILD) && ~isfield(get(CHILD(NumAnnotations+1)),'CData')
%    NumAnnotations=NumAnnotations+1;
%end

if bDrawnow
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bEnoughLoaded, LoadNeedStruct] = EnoughDataLoaded(nSel)
global HISTO;


bEnoughLoaded = 1;
LoadNeedStruct = [];
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadData(bReplot)
%if bResetValue=1, it selects the first item in the box
%otherwise, it tries to keep the same selected item

global HISTO;
MAXGB = 16;

%set(gcf, 'Pointer', 'watch'); drawnow;
try
    if HISTO.bLoadingData~=0
        %another thread is running
        HISTO.bLoadingData = -1-bReplot; %ask it to abort  %-1=don't replot; -2=Replot
        return;
    end
    
    %Set global flag for multi-threading
    HISTO.bLoadingData = 1;
    
    DirNames = {HISTO.SubDir(:).Name};
    
    %set priorities
    nActive = [HISTO.nSel  HISTO.nSel+1  HISTO.nSel+2 HISTO.nSel-1  HISTO.nSel+3:100  HISTO.nSel-2:-1:1];
    nActive(nActive<1 | nActive>length(DirNames)) = [];
    n = nActive(1);  %only load data for the current item.  We can never hope to load all the data even for one of these
    
    if isempty(HISTO.SubDir(n).ChData)
        %FRESH LOADS ARE SPECIAL
        
        [imageData] = bfopenIX83([HISTO.RootDir '\' DirNames{n}], HISTO.Mag);
        if HISTO.bVS200
            imageData = rot90(imageData, -1); %counterclockwize by 90
        end
        absRowY = (1:size(imageData,1))*2^(HISTO.Mag);
        absColX = (1:size(imageData,2))*2^(HISTO.Mag);

        %determine number of channels so future data loads are more efficient
        HISTO.SubDir(n).NumCh = size(imageData,3);
        
        %cache this data (create a new ChData element)
        HISTO.SubDir(n).ChData(1).bPlotted = false;
        HISTO.SubDir(n).ChData(1).mag = HISTO.Mag;
        HISTO.SubDir(n).ChData(1).imageData = imageData;
        HISTO.SubDir(n).ChData(1).absColX = absColX;
        HISTO.SubDir(n).ChData(1).absRowY = absRowY;
        
        %determine default axis
        if 1 %isempty(HISTO.SubDir(n).Axis)
            ImBW = sum(imageData,3);
            ImBW = ImBW-min(ImBW(:));
            %find smallest region that contains data (and eliminate any "gaps")
            RowWithDat = find(sum(ImBW,2)>0); RowWithDat=min(RowWithDat):max(RowWithDat);
            ColWithDat = find(sum(ImBW,1)>0); ColWithDat=min(ColWithDat):max(ColWithDat);
            xlim = absColX(ColWithDat([1 end]));
            ylim = absRowY(RowWithDat([1 end])); %desired x/y lim
            Pos = get(HISTO.UI.Ax, 'Position'); %actual size of axes (in pixels)
            AxW = Pos(3);
            AxH = Pos(4);
            if diff(xlim)/diff(ylim) < AxW/AxH
                %xlim is too narrow -- widen it
                dX = (diff(ylim)*AxW/AxH - diff(xlim))/2;
                %dX = min(dX, xlim(1));  %don't allow negative
                xlim = [xlim(1)-dX  xlim(2)+dX];
                set(HISTO.UI.Ax, 'xlim', xlim);
            elseif diff(ylim)/diff(xlim) < AxH/AxW
                %ylim is too narrow -- widen it
                dY = (diff(xlim)*AxH/AxW - diff(ylim))/2;
                %dY = min(dY, ylim(1));  %don't allow negative
                ylim = [ylim(1)-dY  ylim(2)+dY];
                set(HISTO.UI.Ax, 'ylim', ylim);
            end
            HISTO.SubDir(n).Axis = [];%[xlim ylim]/2^HISTO.Mag;  %save
        end
        
        %redraw
        if bReplot 
            ReplotAx(0,0,1); %note, this has a drawnow. 
            if ThreadAborted
                return;
            end
        end
    end
    
    if ThreadAborted
        return;
    end
    
    %determine total number of slices to analyze
    NumROI = 0;
    for k=1:length(HISTO.SubDir)
        NumROI = NumROI + length(HISTO.SubDir(k).RoiCell);
    end
    UpdateFigName;    
    
    %determine if we need to load any more for this item.
    [bEnoughLoaded, LoadNeedStruct] = EnoughDataLoaded(n);
    while ~bEnoughLoaded
        if ThreadAborted
            return;
        end
        
        %check if we have enough memory for this.
        M=memory;  
        UpdateFigName;
        while length(nActive)>1 && (M.MemUsedMATLAB/1e9 > MAXGB)
            %trying to load higher-priority data, lets clear lower-priority
            
            bClearedSomething = 0;
            while ~bClearedSomething && length(nActive)>1
                if ~isempty(HISTO.SubDir(nActive(end)).ChData)
                    bClearedSomething = 1;
                    HISTO.SubDir(nActive(end)).ChData(:) = []; %clear lowest-priority data
                end
                if isempty(HISTO.SubDir(nActive(end)).ChData)
                    nActive(end) = [];  %clear this
                end
                if ThreadAborted
                    return;
                end
            end
            
            if bClearedSomething
                %recheck memory usage
                M=memory; 
                UpdateFigName;
            end
        end
        
        if (M.MemUsedMATLAB/1e9 > MAXGB)
            %still can't load more data
            waitfor(errordlg('Not enough memory to render current image'));
            ListBoxRefresh;
            HISTO.bLoadingData = 0;
            return;
        end
        
        
        %load the data
        TileX = LoadNeedStruct.NeededTilesXYPix(1,1);
        TileY = LoadNeedStruct.NeededTilesXYPix(1,2);
        disp(LoadNeedStruct.NeededTilesXYPix);
        
        %load data
        [imageData] = bfopenIX83([HISTO.RootDir '\' DirNames{n}], LoadNeedStruct.MaxMag);
        if HISTO.bVS200
            imageData = rot90(imageData, -1); %counterclockwize by 90
        end
        absRowY = (1:size(imageData,1))*2^(HISTO.Mag);
        absColX = (1:size(imageData,2))*2^(HISTO.Mag);
        
        %cache this data (create a new ChData element)
        HISTO.SubDir(n).ChData(end+1).bPlotted = false;
        HISTO.SubDir(n).ChData(end).mag = LoadNeedStruct.MaxMag;
        HISTO.SubDir(n).ChData(end).imageData = imageData;
        HISTO.SubDir(n).ChData(end).absColX = absColX;
        HISTO.SubDir(n).ChData(end).absRowY = absRowY;

        if ThreadAborted
            return;
        end
        
        %redraw
        if bReplot 
            ReplotAx(0,0,1); %note, this has a drawnow. 
            if ThreadAborted
                return;
            end
        end
        
        
        %determine if we need to load any more for this item.
        LoadNeedStruct.NeededTilesXYPix(1,:) = [];
        if isempty(LoadNeedStruct.NeededTilesXYPix)
            bEnoughLoaded = true;
        end
    end
    
    %fill in any missing RoiMean Data
    UpdateRoiAll(n);
    if n==get(HISTO.UI.ListBox,'Value')
        if ThreadAborted
            return;
        end
    end
    
    HISTO.bLoadingData = 0;
catch
    HISTO.bLoadingData = 0;
    rethrow(lasterror);
end

%set(gcf, 'Pointer', 'arrow'); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bAbort = ThreadAborted
drawnow;  %allow UI event to interrupt

if isempty(whos('global', 'HISTO'))
    %program has been exited
    bAbort = 1; %tell the vestigal thread to forget about it
else
    global HISTO
    if HISTO.bLoadingData < 0
        %another thread was called
        bReplot = -HISTO.bLoadingData-1;
        HISTO.bLoadingData = 0;
        LoadData(bReplot); %start over
        bAbort = 1; %tell the vestigal thread to forget about it
    else
        bAbort = 0;  %keep going
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BoxNames = ListBoxRefresh
%if bResetValue=1, it selects the first item in the box
%otherwise, it tries to keep the same selected item

global HISTO;

%Now update the text in the listbox
BoxNames = {HISTO.SubDir(:).Name};
for n=1:length(HISTO.SubDir)
    NumCh = HISTO.SubDir(n).NumCh;
    for ch=1:length(HISTO.SubDir(n).ChData)
        BoxNames{n} = [BoxNames{n}];
    end
end
set(HISTO.UI.ListBox,'String', BoxNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VS200_Callback
global HISTO;

NewVal = get(HISTO.UI.bVS200,'Value');
if HISTO.bVS200 ~= NewVal
    %value changed -- reset all loaded data
    set(gcf, 'Pointer', 'watch'); drawnow;
    HISTO.bVS200 = NewVal;
    
    for n=1:length(HISTO.SubDir)
        HISTO.SubDir(n).ChData = [];
    end
    ListBox_Callback;
    
    set(gcf, 'Pointer', 'arrow'); drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HISTO.SubDir(n).bFlipLR = get(HISTO.UI.bFlipLR, 'value');
if HISTO.SubDir(n).bFlipLR
    %'reverse'
    set(HISTO.UI.Ax,'xdir','reverse');
else
    set(HISTO.UI.Ax,'xdir','normal');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FlipLR_Callback
global HISTO;

n = HISTO.nSel;
try
    HISTO.SubDir(n).bFlipLR = get(HISTO.UI.bFlipLR, 'value');
    if HISTO.SubDir(n).bFlipLR
        %'reverse'
        set(HISTO.UI.Ax,'xdir','reverse');
    else
        set(HISTO.UI.Ax,'xdir','normal');
    end
catch
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mag_Callback
global HISTO;

NewMag = get(HISTO.UI.Mag,'Value');
if HISTO.Mag ~= NewMag
    %Mag changed -- reset all loaded data
    set(gcf, 'Pointer', 'watch'); drawnow;
    HISTO.Mag = NewMag;
    
    for n=1:length(HISTO.SubDir)
        HISTO.SubDir(n).ChData = [];
    end
    OldAxLimits = axis(HISTO.UI.Ax);
    ListBox_Callback;
    axis(HISTO.UI.Ax, OldAxLimits);
    
    set(gcf, 'Pointer', 'arrow'); drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ListBox_Callback
global HISTO;
set(gcf, 'Pointer', 'watch'); drawnow;

MainFigure_ResizeFcn;
if ~HISTO.bFirstClick
    %save axis/zoom info for previously selected item
    axes(HISTO.UI.Ax);
    HISTO.SubDir(HISTO.nSel).Axis = []; %round(axis(HISTO.UI.Ax));
end
HISTO.bFirstClick = false;

%set new selected item
HISTO.nSel = get(HISTO.UI.ListBox,'Value');

%clear flags -- all images must be redrawn
for n=1:length(HISTO.SubDir)
    for k=1:length(HISTO.SubDir(n).ChData)
        HISTO.SubDir(n).ChData(k).bPlotted = false;
    end
    cla(HISTO.UI.Ax);  %clear.  NOTE: this needs fix for ROI
end

%set(gcf, 'Pointer', 'watch'); drawnow; 
%set(gcf, 'Pointer', 'arrow'); drawnow;
bLoadDataCalled = ReplotAx(1,1,1);
if ~bLoadDataCalled
    LoadData(0);
end

set(gcf, 'Pointer', 'arrow'); drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RGBCheck_Callback

Edit_Callback;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SingleRGB(RGBhit)
global HISTO

Flds = {'R' 'G' 'B'};
Hit = [0 0 0];
Check = [0 0 0];
%first asses the situation
for c=1:3
    Check(c) = get(HISTO.UI.RGBCheck.(Flds{c}), 'Value');
    if Flds{c}==RGBhit
        Hit(c) = 1;
    end
end
if all(Hit==Check)
    %untoggle -- turn everything on
    VAL = [1 1 1];
else
    %turn only this color on
    VAL = Hit;
end
for c=1:3
    set(HISTO.UI.RGBCheck.(Flds{c}), 'Value', VAL(c));
end
RGBCheck_Callback;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Edit_Callback
global HISTO

set(gcf, 'Pointer', 'watch'); drawnow;

OldCh = [HISTO.DispRGB(:).Ch];
OldAxLimits = axis(HISTO.UI.Ax);

Flds = {'Ch' 'Min' 'Max' 'Gam'};
for k=1:length(Flds)
    for c=1:3
        tmp = get(HISTO.UI.RGB(c).(Flds{k}), 'string');
        num = str2num(tmp);
        if ~isempty(num)
            tmp = num;
        end
        HISTO.DispRGB(c).(Flds{k}) = tmp;
    end
end
Flds = {'R' 'G' 'B'};
for c=1:3
    HISTO.DispRGB(c).On = get(HISTO.UI.RGBCheck.(Flds{c}), 'Value');
end
%all images must be redrawn
for n=1:length(HISTO.SubDir)
    for k=1:length(HISTO.SubDir(n).ChData)
        HISTO.SubDir(n).ChData(k).bPlotted = false;
    end
    %cla(HISTO.UI.Ax);  %clear.  NOTE: this needs fix for ROI
end
%set(gcf, 'Pointer', 'watch'); drawnow;

%MainFigure_ResizeFcn;
%axes(HISTO.UI.Ax); 
%HISTO.SubDir(HISTO.nSel).Axis = []; %round(axis(HISTO.UI.Ax));
ReplotAx(1,0,0);
axis(HISTO.UI.Ax, OldAxLimits);

uicontrol(HISTO.UI.ListBox); %set focus to ListBox
set(gcf, 'Pointer', 'arrow'); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonREDRAW

global HISTO;

set(gcf, 'Pointer', 'watch'); drawnow;
MainFigure_ResizeFcn;
axes(HISTO.UI.Ax); 
HISTO.SubDir(HISTO.nSel).Axis = []; %round(axis(HISTO.UI.Ax));
OldAxLimits = axis(HISTO.UI.Ax);
ReplotAx(0,1,1);
axis(HISTO.UI.Ax, OldAxLimits);
set(gcf, 'Pointer', 'arrow'); drawnow;

UpdateFigName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UpdateFigName
global HISTO;
%determine total number of slices to analyze
NumROI = 0;
for k=1:length(HISTO.SubDir)
    NumROI = NumROI + length(HISTO.SubDir(k).RoiCell);
end
M=memory;  set(gcf, 'Name', ['HistoIX83 ' num2str(M.MemUsedMATLAB/1e9) ' GB.   #ROI = ' num2str(NumROI)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonNewCell_Callback(TypeNuc, TypeCyto)
global HISTO

%if isempty(HISTO.SubDir(HISTO.nSel).RoiBkg)
%    waitfor(errordlg('Please draw background ROI before any Cell ROIs'));
%    return;
%end

try
    jCell = length(HISTO.SubDir(HISTO.nSel).RoiCell) + 1;
    NewCell.bShow = true;
    set(gcf, 'Pointer', 'watch'); drawnow;  %will be replaced by pan function
    NewCell.Nuc  = NewROI(TypeNuc,  'm', HISTO.nSel, jCell, 'Nuc');
    set(gcf, 'Pointer', 'watch'); drawnow;  %will be replaced by pan function
    NewCell.Cyto = NewROI(TypeCyto, 'g', HISTO.nSel, jCell, 'Cyto');
    set(gcf, 'Pointer', 'arrow');
    %pan(gca,'on')  %hand
    
    HISTO.SubDir(HISTO.nSel).RoiCell = [HISTO.SubDir(HISTO.nSel).RoiCell  NewCell];
    uicontrol(HISTO.UI.ListBox); %set focus to ListBox
    
    %calc Means
    UpdateRoi(HISTO.nSel, jCell, 'Nuc', []);
    UpdateRoi(HISTO.nSel, jCell, 'Cyto', []);
    
    UpdateFigName;
catch
    waitfor(errordlg(lasterr));
    set(gcf, 'Pointer', 'arrow');
    uicontrol(HISTO.UI.ListBox); %set focus to ListBox
    ListBox_Callback;  %redraw
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MouseButtonDown
global HISTO
%if isstruct(HISTO.UI.RoiMovement) || HISTO.UI.RoiMovement~=0
%    warning('MouseButtonDown occured with unusual UI.RoiMovement Value');
%end
HISTO.UI.RoiMovement = 1;  %allow new movement.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RoiMovementUpdate(nDir, jCell, NucCyto, Pos)
%keep track of position until done editing
global HISTO
if isstruct(HISTO.UI.RoiMovement) || HISTO.UI.RoiMovement==1
    %allow new movement update
    HISTO.UI.RoiMovement = [];
    HISTO.UI.RoiMovement.nDir = nDir;
    HISTO.UI.RoiMovement.jCell = jCell;
    HISTO.UI.RoiMovement.NucCyto = NucCyto;
    HISTO.UI.RoiMovement.Pos = Pos;
else  %HISTO.UI.RoiMovement==0
    %movement after a button up command -- reset the pos to the previous one
    if jCell==0
        ROI = HISTO.SubDir(nDir).RoiBkg;
    else
        ROI = HISTO.SubDir(nDir).RoiCell(jCell).(NucCyto);
    end
    %kill this ROI, and create a new one with the correct position.
    ROI = DestructROI(ROI);
    ROI = ReloadROI(ROI, nDir, jCell, NucCyto);
    if jCell==0
        HISTO.SubDir(nDir).RoiBkg = ROI;
    else
        HISTO.SubDir(nDir).RoiCell(jCell).(NucCyto) = ROI;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MouseButtonUp
global HISTO
if isstruct(HISTO.UI.RoiMovement)
    M = HISTO.UI.RoiMovement;
    UpdateRoi(M.nDir, M.jCell, M.NucCyto, M.Pos);
    HISTO.UI.RoiMovement = 0;    %done editing
    uicontrol(HISTO.UI.ListBox); %set focus to ListBox
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateRoiAll(nDir)
global HISTO
if ~isempty(HISTO.SubDir(nDir).RoiBkg)
    UpdateRoi(nDir, 0, 0, []);
end
for j = 1:length(HISTO.SubDir(nDir).RoiCell)
    UpdateRoi(nDir, j, 'Nuc', []);
    %UpdateRoi(nDir, j, 'Cyto', []);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateRoi(nDir, jCell, NucCyto, Pos, bForce)

global HISTO

if jCell==0
    ROI = HISTO.SubDir(nDir).RoiBkg;
else
    ROI = HISTO.SubDir(nDir).RoiCell(jCell).(NucCyto);
end
if isempty(ROI)
    return;
end

try
    bForce;
catch
    bForce = 0;
end

%store new position
if ~isempty(Pos)
    ROI.Pos = Pos;
end

%%%CODE ADDED 3_24_15 --%NO NEED TO CALCULATE IN REAL TIME -- THIS PROGRAM IS MEANT FOR A SLOW OFF-LINE CALCULATION
if jCell==0
    HISTO.SubDir(nDir).RoiBkg = ROI;
else
    HISTO.SubDir(nDir).RoiCell(jCell).(NucCyto) = ROI;
end
return;  
%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if bForce || isempty(ROI.Mean) ||  ~isempty(Pos)
    %clear the Mean data
    ROI.Mean = nan(1, HISTO.SubDir(nDir).NumCh);
end

%calculate Means only if images are loaded, and some mean hasn't been calculated yet
if ~isempty(HISTO.SubDir(nDir).ChData) && ~isempty(HISTO.SubDir(nDir).ChData(1).Im) &&  any(isnan(ROI.Mean(:)))
    
    %generate mask
    if strcmpi(ROI.Type,'imblob')
        BW = false(size(HISTO.SubDir(nDir).ChData(1).Im));
        BW(ROI.Pix) = true;
    elseif ~isempty(ROI.h)
        BW = createMask(ROI.h);
    elseif strcmpi(ROI.Type, 'imellipse')
        %get equation parameters for ellipse
        Pos = ROI.Pos;
        x0 = Pos(1) + Pos(3)/2;
        y0 = Pos(2) + Pos(4)/2;
        a = Pos(3)/2;
        b = Pos(4)/2;
        %convert this to a mask
        [xx,yy] = meshgrid(1:size(HISTO.SubDir(nDir).ChData(1).Im,2),1:size(HISTO.SubDir(nDir).ChData(1).Im,1));
        BW = ((xx-x0).^2/a^2 + (yy-y0).^2/b^2 < 1); %points within the ellipse
        
        %draw a quick ellipse
        h = ellipse(a,b,0,x0,y0,'w');
        set(h,'Parent', HISTO.UI.Ax)
        drawnow;
    else
        %not the image being viewed.
        tmpH = feval(ROI.Type, HISTO.UI.Ax, ROI.Pos);
        BW = createMask(tmpH);  %warning: this may cause problems if image sizes change
        delete(tmpH);
    end
    
    %calculate means
    for ch=1:HISTO.SubDir(nDir).NumCh
        try
            if isnan(ROI.Mean(ch))
                tmp = HISTO.SubDir(nDir).ChData(ch).Im;
                ROI.Mean(z,ch) = mean(tmp(BW));
            end
        catch
        end
    end
    if jCell==0
        HISTO.SubDir(nDir).RoiBkg = ROI;
    else
        HISTO.SubDir(nDir).RoiCell(jCell).(NucCyto) = ROI;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ROI = NewROI(Type, color, nDir, jCell, NucCyto)
global HISTO

ROI.Type = Type;
ROI.h = feval(ROI.Type, HISTO.UI.Ax);
ROI.h.setColor(color);
ROI.id = addNewPositionCallback(ROI.h, @(pos) RoiMovementUpdate(nDir,jCell,NucCyto,pos));
ROI.Pos = ROI.h.getPosition;
ROI.Mean = [];

BW = createMask(ROI.h);
if sum(BW(:)) == 0
    error('ROI has zero area');
end

%set uicontextmenu
ui=uicontextmenu;
if jCell>0
    uimenu(ui,'label','Delete','callback',CallbackString('HistoIX83','''ROICheck_ButtonDown''', jCell));
end
set(ROI.h, 'UIContextMenu', ui);
set(get(ROI.h,'children'), 'UIContextMenu', ui);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ROI = ReloadROI(ROI, nDir, jCell, NucCyto)
global HISTO

if jCell==0
    color = 'b';
elseif strcmpi(NucCyto,'Nuc')
    color = 'm';
else
    color = 'g';
end

if strcmpi(ROI.Type,'imblob')
    ROI.h = line(ROI.Pos([1:end 1],1),ROI.Pos([1:end 1],2));
    set(ROI.h,'color','w', 'linewidth', 2, 'parent', HISTO.UI.Ax);
    ROI.id = [];
else
    ROI.h = feval(ROI.Type, HISTO.UI.Ax, ROI.Pos);
    ROI.h.setColor(color);
    ROI.id = addNewPositionCallback(ROI.h, @(pos) RoiMovementUpdate(nDir,jCell,NucCyto,pos));
end

%set uicontextmenu
ui=uicontextmenu;
if jCell>0
    uimenu(ui,'label','Delete','callback',CallbackString('HistoIX83','''ROICheck_ButtonDown''', jCell));
end
set(ROI.h, 'UIContextMenu', ui);
set(get(ROI.h,'children'), 'UIContextMenu', ui);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReplotROI(bDrawnow)
global HISTO;

%regenerate all ROI's
%if ~isempty(HISTO.SubDir(HISTO.nSel).RoiBkg)
%    HISTO.SubDir(HISTO.nSel).RoiBkg = ReloadROI(HISTO.SubDir(HISTO.nSel).RoiBkg, HISTO.nSel, 0, 0);
%end
for j=1:length(HISTO.SubDir(HISTO.nSel).RoiCell)
    if HISTO.SubDir(HISTO.nSel).RoiCell(j).bShow
        HISTO.SubDir(HISTO.nSel).RoiCell(j).Nuc  = ReloadROI(HISTO.SubDir(HISTO.nSel).RoiCell(j).Nuc,  HISTO.nSel, j, 'Nuc');
        HISTO.SubDir(HISTO.nSel).RoiCell(j).Cyto = ReloadROI(HISTO.SubDir(HISTO.nSel).RoiCell(j).Cyto, HISTO.nSel, j, 'Cyto');
    end
end
if bDrawnow
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ClearROIHandles
global HISTO;

for n=1:length(HISTO.SubDir)
    try
        HISTO.SubDir(n).RoiBkg.Pos = HISTO.SubDir(n).RoiBkg.h.getPosition;
    catch
    end
    HISTO.SubDir(n).RoiBkg = DestructROI(HISTO.SubDir(n).RoiBkg);
    
    for j=1:length(HISTO.SubDir(n).RoiCell)
        try
            HISTO.SubDir(n).RoiCell(j).Nuc.Pos = HISTO.SubDir(n).RoiCell(j).Nuc.h.getPosition;
        catch
        end
        HISTO.SubDir(n).RoiCell(j).Nuc = DestructROI(HISTO.SubDir(n).RoiCell(j).Nuc);
        
        try
            HISTO.SubDir(n).RoiCell(j).Cyto.Pos = HISTO.SubDir(n).RoiCell(j).Cyto.h.getPosition;
        catch
        end
        HISTO.SubDir(n).RoiCell(j).Cyto = DestructROI(HISTO.SubDir(n).RoiCell(j).Cyto);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ROI = DestructROI(ROI)
if ~isempty(ROI) && isstruct(ROI)
    warning off;
    if ~isempty(ROI.h) && ~isempty(ROI.id)
        try
            removeNewPositionCallback(ROI.h, ROI.id);
        catch
        end
    end
    if ~isempty(ROI.h)
        try
            delete(ROI.h);
        catch
        end
    end
    ROI.id = [];
    ROI.h = [];
    warning on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ROI = UnlinkROI(ROI)
warning off;
ROI.id = [];
ROI.h = [];
warning on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteAllRoi
%delte all Roi for current selected image
global HISTO
%right-click = delete ROI
n = HISTO.nSel;
for j=length(HISTO.SubDir(n).RoiCell):-1:1
    try
        delete(HISTO.SubDir(n).RoiCell(j).Nuc.h);
    catch
    end
    try
        delete(HISTO.SubDir(n).RoiCell(j).Cyto.h);
    catch
    end
    HISTO.SubDir(n).RoiCell(j) = [];
end
ClearROIHandles; %clear ROI UI objects
ReplotROI(1);   %replot ROI objects
uicontrol(HISTO.UI.ListBox); %set focus to ListBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ROICheck_ButtonDown(j)
global HISTO
%right-click = delete ROI
n = HISTO.nSel;
if j<=length(HISTO.SubDir(n).RoiCell)
    try
        delete(HISTO.SubDir(n).RoiCell(j).Nuc.h);
    catch
    end
    try
        delete(HISTO.SubDir(n).RoiCell(j).Cyto.h);
    catch
    end
    HISTO.SubDir(n).RoiCell(j) = [];
    
    ClearROIHandles; %clear ROI UI objects
    ReplotROI(1);   %replot ROI objects
end
uicontrol(HISTO.UI.ListBox); %set focus to ListBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RGB] = MakeRGB(imageData, DispRGB)
%DispRGB
%  .Ch  - defines which is R, G, B
%  .Min,Max,Gam =  [min max gamma]   note; if max<=1, then image will first be normalized [0 1].

tmpL = imageData(:,1:round(end/2),1);
tmpR = imageData(:,round(end/2)+1:end,1);

RGB = zeros(size(imageData,1), size(imageData,2), 3, 'uint8'); %just for display, so low-order dat is fine
for k=1:3 %RGB
    try
        if DispRGB(k).On
            tmp = imageData(:,:,DispRGB(k).Ch);
        else
            tmp = imageData(:,:,1)*0;
        end
    catch
        tmp = imageData(:,:,1)*0;
    end
    RGB(:,:,k) = uint8(255*Normalize(tmp, [DispRGB(k).Min  DispRGB(k).Max  DispRGB(k).Gam]));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Im = Normalize(Im, ScaleInfo)
%ScaleInfo = [min max gamma]
%note; if max<=1, then image will first be normalized [0 1].

Im = double(Im);
if ScaleInfo(2)<=1
    Im = Im - min(Im(:));
    Im = Im/max(Im(:));
end
%if ScaleInfo(1)==0
%    Im = Im - min(Im(:));
%end
%normalize with given min/max
Im = (Im - ScaleInfo(1))/(ScaleInfo(2)-ScaleInfo(1));
%clip
Im(Im<0) = 0;
Im(Im>1) = 1;

%now apply gamma
Im = Im.^ScaleInfo(3);

%Im = imtophat(Im,strel('disk',100));
%Im = Im/max(Im(:));

return;
se = strel('ball',10,1);
Im = imsubtract(imadd(Im,imtophat(Im,se)), imbothat(Im,se));
disp([min(Im(:)) max(Im(:))]);
Im = Im - min(Im(:));
Im = Im/max(Im(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = CallbackString(RootFunction, varargin)
%function str = CallbackString(RootFunction, varargin)
%
%generates a string to use as a callback funciton
%the arguments can be function names, strings, or numbers
%NOTE: function names are in single quotes, strings in triple quotes, or in cell-brackets
%Example: CallbackString(MyFunction, 'gcbo', '''hello''', {'bracketed'}, 3) will return (as a string:)
%        MyFunction(gcbo,'hello','bracketed',3)

str = strcat(RootFunction, '(');
for k=1:nargin-1
    if iscell(varargin{k})
        str = strcat(str, '''', char(varargin{k}), '''');    %we must add the quotes around this string
    elseif ischar(varargin{k})
        str = strcat(str, varargin{k});
    else
        %this is numerical, or matrix input -- must convert it into a string.
        %this requires special attention if there are multiple rows
        A = varargin{k};
        if isempty(A)
            str = strcat(str,'[]');
        else
            str = strcat(str, '[');  %opening bracket
            str = strcat(str, num2str(A(1,:)));  %first row
            for nRow=2:size(A,1)
                str = strcat(str, ';', num2str(A(nRow,:)));  %subsequent rows
            end
            str = strcat(str, ']');  %closing bracket
        end
    end
    if k<nargin-1
        str = strcat(str, ',');  %comma
    end
end
str = strcat(str, ')');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


