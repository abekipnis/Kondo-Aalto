function Header = ImportVERTHeader(filename)
%IMPORTVERTHEADER Import header data from a Createc VERT file
% Read out all header information until the max header byte of 16384 bytes.
% Input: full path to the VERT data file
% Output: struct with data fields and their values
% Based on code written by Benjamin Alldritt

%% Open the file
fid = fopen(filename);
if fid == -1
    error('File not found or can not be read!')
end

%% Check VERT file version
filetype = fgetl(fid);
if strcmp('[ParVERT32]',filetype)
    STMAFMVersion = 4;
    headerBytes = 16372;
elseif strcmp('[ParVERT30]',filetype)
    STMAFMVersion = 3;
    headerBytes = 16371;
% elseif strcmp('[Paramco32]',filetype)
%     % Both versions 3 and 4 have this first row. Read more.
%     fgetl(fid); ft34 = fgetl(fid);
%     if contains(ft34, 'Titel(File)')
%         STMAFMVersion = 4;
%     else
%         STMAFMVersion = 3;
%     end
%     headerBytes = 16384 - 13 - 16 - length(ft34) - 2; 
else
% If none of these is found, stop and signal file error.
    error('Createc VERT file version does not match')
end


%% Start reading the header
header = fread(fid,[1,headerBytes],'*char');

% First, insert a carriage return to get rid of incorrectly formatted dates at end of file,
% split the header into separate lines, then remove garbage filler data at
% the end of the header. Strip the leading and tailing whitespace.
% Finally split the values and remove leading and tailing whitespace.
SplittedLine = insertBefore(header, 'DSP-COMPDATE', char(13));
SplittedLine = string(splitlines(SplittedLine));
SplittedLine = SplittedLine(1:find(contains(SplittedLine,'DSP-COMPDATE')),:);
%SplittedLine = strip(SplittedLine);

% Knife for version 4: one line is not delimited by " = "
if STMAFMVersion == 4
    SplittedLine(cellfun('isempty', SplittedLine)) = [];
end

% Sometimes the memo does not contain a "=": take the content and append it
% to the previous line
noDelim = SplittedLine(~contains(SplittedLine, '='));
if ~isempty(noDelim)
    memoIndex = find(SplittedLine == noDelim) -1;
    SplittedLine(memoIndex) = SplittedLine(memoIndex) + ' ' + noDelim;
    SplittedLine(memoIndex +1) = [];
end

SplittedLine = split(SplittedLine, '=');
SplittedLine = strip(SplittedLine);

% Header Values
% Read in header values from the split data.
Header.ScanPixels_X = str2double(SplittedLine(find(contains(SplittedLine,'Num.X / Num.X')),2));
Header.ScanPixels_Y = str2double(SplittedLine(find(contains(SplittedLine,'Num.Y / Num.Y')),2));
Header.GainX = str2double(SplittedLine(find(contains(SplittedLine,'GainX / GainX')),2));
Header.GainY = str2double(SplittedLine(find(contains(SplittedLine,'GainY / GainY')),2));
Header.GainZ = str2double(SplittedLine(find(contains(SplittedLine,'GainZ / GainZ')),2));
Header.GainPreamplifier = str2double(SplittedLine(find(contains(SplittedLine,'Gainpreamp / GainPre 10^')),2));
Header.ChannelCount = str2double(SplittedLine(find(contains(SplittedLine,'Channels / Channels')),2));
Header.VertManGain = str2double(SplittedLine(find(contains(SplittedLine,'Vertmangain')),2));

% These seem to be useless nonsense
Header.DACToZConversionFactor = str2double(SplittedLine(find(contains(SplittedLine,'Dacto[A]z')),2));
Header.DACToXYConversionFactor = str2double(SplittedLine(find(contains(SplittedLine,'Dacto[A]xy')),2));

Header.ScanRange_X = str2double(SplittedLine(find(contains(SplittedLine,'Length x[A]')),2))*10^-10;
Header.ScanRange_Y = str2double(SplittedLine(find(contains(SplittedLine,'Length y[A]')),2))*10^-10;
Header.ScanOffset_X = str2double(SplittedLine(find(contains(SplittedLine,'Scanrotoffx / OffsetX')),2));
Header.ScanOffset_Y = str2double(SplittedLine(find(contains(SplittedLine,'Scanrotoffy / OffsetY')),2));
Header.Bias = str2double(SplittedLine(find(contains(SplittedLine,'Biasvolt[mV]')),2));
Header.Current = str2double(SplittedLine(find(contains(SplittedLine,'Current[A]')),2));

Header.PiezoX = str2double(SplittedLine(find(contains(SplittedLine,'Xpiezoconst / Xpiezoconst')),2));
Header.PiezoY = str2double(SplittedLine(find(contains(SplittedLine,'YPiezoconst / YPiezoconst')),2));
Header.PiezoZ = str2double(SplittedLine(find(contains(SplittedLine,'ZPiezoconst / ZPiezoconst')),2));

Header.LockinFreq = str2double(SplittedLine(find(contains(SplittedLine,'LockinFreq')),2));
Header.LockinRC = str2double(SplittedLine(find(contains(SplittedLine,'LockinRC')),2));
Header.LockinPhase = str2double(SplittedLine(find(contains(SplittedLine,'LockinPhase')),2));
Header.LockinModulation = str2double(SplittedLine(find(contains(SplittedLine,'LockinAmpl')),2));
Header.ZOffset_DAC = str2double(SplittedLine(find(contains(SplittedLine,'Zoffset')),2));
Header.ZOffset_A = Header.ZOffset_DAC*Header.PiezoZ*1e-3;

% Latest saved scan file before the spectrum:
Header.LinkFileName = convertStringsToChars(SplittedLine(find(contains(SplittedLine,'LinkFilename / LinkFilename')),2));

% These guys are saved later
% Header.VertSpecPosX = str2double(SplittedLine(find(contains(SplittedLine,'VertSpecPosX')),2));
% Header.VertSpecPosY = str2double(SplittedLine(find(contains(SplittedLine,'VertSpecPosY')),2));

% V4 knife: Current is saved under a different name
if STMAFMVersion == 4
    Header.Current = str2double(SplittedLine(find(contains(SplittedLine,'FBLogIset')),2));
    % Also save the bwd/fwd setting
    Header.BackwardForward = str2double(SplittedLine(find(contains(SplittedLine,'Scanmode / ScanXMode')),2));
end

Header.ACQ_Time = str2double(SplittedLine(find(contains(SplittedLine,'Sec/Image:')),2));
Header.ScanAngle = str2double(SplittedLine(find(contains(SplittedLine,'Rotation / Rotation')),2));
Header.ZControllerSetpoint = str2double(SplittedLine(find(contains(SplittedLine,'FBLogIset')),2));
Header.ZControllerIntegralGain = str2double(SplittedLine(find(contains(SplittedLine,'FBIntegral')),2));
Header.ZControllerProportionalGain = str2double(SplittedLine(find(contains(SplittedLine,'FBProp')),2));

Header.CHMode = str2double(SplittedLine(find(contains(SplittedLine,'CHMode / CHMode')),2));

Header.LHeTemperature = str2double(SplittedLine(find(contains(SplittedLine,'T_AUXADC7[K]')),2));
Header.STMTemperature = str2double(SplittedLine(find(contains(SplittedLine,'T_AUXADC6[K]')),2));

%% Then beyond the header

for jj = 1:1
    SkipLine = fgetl(fid);
end

if STMAFMVersion == 4
    % Read the mystery line
    InfoLine = fgetl(fid);
    
    % Split the mystery line into nice cell array
    Info = split(InfoLine);
    
    % Save parameters to Header
    Header.NPoints = str2double(Info{2});
    Header.VertPosX_DAC = str2double(Info{3}); % Apparently in DAC relative to scan offset
    Header.VertPosY_DAC = str2double(Info{4});
    Header.ChannelList = str2double(Info{5}); % Sum of binary channel values, see help file
    Header.ChannelCount = str2double(Info{6}); % ???
    Header.OutChannelList = str2double(Info{7}); % ???
    Header.OutChannelCount = str2double(Info{8}); % Determines which columns are which
    Header.XPos_nm = str2double(Info{9}); % Spectrum X-coordinate in the approach area
    Header.YPos_nm = str2double(Info{10}); % Spectrum Y-coordinate in the approach area
end

fclose(fid);

end