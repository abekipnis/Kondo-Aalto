function [Header, output, Im] = ReadDAT(file)
% Read Createc DAT files, irrespective of software version
% Written by Ben Alldritt, patched for Createc v4 by Markus Aapro
% Inputs:
%   - file: a full file path to the .dat file
% Outputs:
%   - Header: a struct containing field values from the header
%   - output: the scan data channels in a cell array: first half are forward
%   - Im: output-array converted into an SPIW-compatible struct

% Open the Createc DAT file. If nothing, return an error.
fid = fopen(file);
if fid == -1
    error('File not found or can not be read!')
end

% DAT files contain in first row "[Paramet32]" in Version 2 and
% "[Parameter]" in older Versions. The newest Version 3 can
% compress the files, which is marked by "[Paramco32]". 13 bytes long
filetype = fgetl(fid);

if strcmp('[Parameter]',filetype)

    STMAFMVersion = 1;
    headerBytes = 16371;

elseif strcmp('[Paramet32]',filetype)

    STMAFMVersion = 2;
    headerBytes = 16371;

elseif strcmp('[Paramco32]',filetype)

    % Both versions 3 and 4 have this first row. Read another line.
    fgetl(fid); ft34 = fgetl(fid);
    if contains(ft34, 'Titel(File)')
        STMAFMVersion = 4;
    else
        STMAFMVersion = 3;
    end
    headerBytes = 16384 - 13 - 16 - length(ft34) - 2; 

else

    % If none of these is found, stop and signal file error.
    error('Createc DAT file version does not match')

end

% Read out all header information until the max header byte of 16384 bytes.
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
Header.DACToZConversionFactor = str2double(SplittedLine(find(contains(SplittedLine,'Dacto[A]z')),2));
Header.DACToXYConversionFactor = str2double(SplittedLine(find(contains(SplittedLine,'Dacto[A]xy')),2));
Header.ScanRange_X = str2double(SplittedLine(find(contains(SplittedLine,'Length x[A]')),2))*10^-10;
Header.ScanRange_Y = str2double(SplittedLine(find(contains(SplittedLine,'Length y[A]')),2))*10^-10;
Header.ScanOffset_X = str2double(SplittedLine(find(contains(SplittedLine,'Scanrotoffx / OffsetX')),2));
Header.ScanOffset_Y = str2double(SplittedLine(find(contains(SplittedLine,'Scanrotoffy / OffsetY')),2));
Header.Bias = str2double(SplittedLine(find(contains(SplittedLine,'Biasvolt[mV]')),2));
Header.Current = str2double(SplittedLine(find(contains(SplittedLine,'Current[A]')),2));

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
Header.PiezoX = str2double(SplittedLine(find(contains(SplittedLine,'Xpiezoconst')),2));
Header.CHMode = str2double(SplittedLine(find(contains(SplittedLine,'CHMode / CHMode')),2));

if STMAFMVersion == 1
    BytePerPixel = 2;
    % Header + 2 unused "NULL"-Bytes
    fseek(fid,2,'cof');
elseif STMAFMVersion == 2
    BytePerPixel = 4;
    % Header + 4 unused "NULL"-Bytes
    fseek(fid,4,'cof');
elseif STMAFMVersion == 3
    BytePerPixel = 4;
    % No Seek of additional bytes, since they are compressed:
elseif STMAFMVersion == 4
    BytePerPixel = 4;
    % No idea about additional bytes: see the Createc software for file
    % structure
else
    % Nothing -> Case is already catched by returning an empty image.
end

% Continue from end of header position, then readout the binary data. In this
% case, the binary data is compressed (zlib compression method).
data = fread(fid);
fclose(fid);

% After reading the raw binary data, convert it to uint8, then decode the
% entire raw data.
output = uint8(data);
output = transpose(zlibdecode(output));

% Convert to 4 byte value (since we are using compressed data for this
% test.
output = typecast(output, 'uint32');

% Remove starting 4 byte buffer
output(1,:) = [];

% Convert to floating point values
output = typecast(output, 'single');

% Convert column vector matrix ScanPixels_X wide
output = vec2mat(output,Header.ScanPixels_X);

% Remove last 4 rows (buffer)
if size(output,1) > 1 
    output(length(output)-3:end,:) = []; 
end

% Knife for version 4:
if STMAFMVersion ==4
    output(end,:) = [];
end

% Convert to cell array based on channel number.

% %Knife for V4: separate bwd/fwd channels
% if STMAFMVersion == 4
% %     output = mat2cell(output, ...
% %         Header.ScanPixels_Y*ones(1,Header.BackwardForward*Header.ChannelCount), ...
% %         Header.ScanPixels_X);
% %     % This tries to create a [dataChannels, direction]-array, but
% %     %  ChannelCount is actually dataChannels*directions
% %     output = mat2cell(output, ...
% %     Header.ScanPixels_Y*ones(1,Header.ChannelCount), ...
% %     [Header.ScanPixels_X, Header.BackwardForward]);
% else
%     output = mat2cell(output, ...
%         Header.ScanPixels_Y*ones(1,Header.ChannelCount), ...
%         Header.ScanPixels_X);
% end
output = mat2cell(output, ...
        Header.ScanPixels_Y*ones(1,Header.ChannelCount), ...
        Header.ScanPixels_X);

% Test converting to SPIW data format (not used in Aalto Atoms code)
Im.data = output;
Im.points = Header.ScanPixels_X;
Im.lines = Header.ScanPixels_Y;
Im.width = Header.ScanRange_X;
Im.height = Header.ScanRange_Y;

end


function output = zlibdecode(input)
%ZLIBDECODE Decompress input bytes using ZLIB.
%
%    output = zlibdecode(input)
%
% The function takes a compressed byte array INPUT and returns inflated
% bytes OUTPUT. The INPUT is a result of GZIPENCODE function. The OUTPUT
% is always an 1-by-N uint8 array. JAVA must be enabled to use the function.
%
% See also zlibencode typecast

error(nargchk(1, 1, nargin));
error(javachk('jvm'));
if ischar(input)
  warning('zlibdecode:inputTypeMismatch', ...
          'Input is char, but treated as uint8.');
  input = uint8(input);
end
if ~isa(input, 'int8') && ~isa(input, 'uint8')
    error('Input must be either int8 or uint8.');
end

buffer = java.io.ByteArrayOutputStream();
zlib = java.util.zip.InflaterOutputStream(buffer);
zlib.write(input, 0, numel(input));
zlib.close();
output = typecast(buffer.toByteArray(), 'uint8')';


end