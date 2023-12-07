function [header, V, Z, fullData, channelData] = ImportSpecGrid(options)
    arguments
        options.file = 'no';
        options.path = 'no';
        options.channel = 2;
        options.normalisation string = 'median';
    end
%IMPORTSPECGRID Read, normalise and import specgrid data
%   Optional inputs:
%       - file: .specgrid file to import
%       - path: path to .specgrid file. Launches file browser by default.
%       - channel: integer specifying the data channel to be imported into
%       channelData variable (default 2 for dI/dV)
%       - normalisation: how to normalise the channelData
%           - raw (default): no normalisation
%           - median: divide by pointwise median
%           - Imean: divide by pointwise average current
%           - IV: divide by the pointwise IV vector
% Based on code supplied by Gerhard Meyer and Viliam Vano
% FIXME: broken dx, dy
% FIXME: broken VZ table

% File base and path
if strcmp(options.file, 'no')
    [file, path] = uigetfile('*.specgrid',...
        'Select the specgrid to import', ...
        'Multiselect', 'off');
else
    file = options.file;
    path = options.path;
end

% Exit in case no files were imported
if isfloat(file)
    disp('No file selected: exiting')
    return;
end

%% Read and import

% Twiddle with float formats
formatOrig = format;
if ~strcmp(formatOrig.NumericFormat, 'longE')
    format longE
end

filename = fopen([path file]);
if filename == -1
    error('File not found or can not be read!')
end

header.ver =                        fread(filename, 1,  'int32');
header.nx =                         fread(filename, 1,  'int32');
header.ny =                         fread(filename, 1,  'int32');
% FIXME: broken hardcoded piezo constants
header.dx =                         fread(filename, 1,  'int32')*83.5/38;%*10V/524287*gain*xy piezo constant
header.dy =                         fread(filename, 1,  'int32')*83.5/38;

header.specxgrid =                  fread(filename, 1,  'int32'); 
header.specygrid =                  fread(filename, 1,  'int32');
header.vertpoints =                 fread(filename, 1,  'int32');
header.vertmandelay =               fread(filename, 1,  'int32');
header.vertmangain =                fread(filename, 1,  'int32');
header.biasvoltage =                fread(filename, 1,  'single');
header.tunnelcurrent =              fread(filename, 1,  'single');
header.imagedatasize =              fread(filename, 1,  'int32'); % Almost the file size of the .specgrid.dat
header.specgriddatasize =           fread(filename, 1,  'int32'); % Almost the file size of the .specgrid
header.specgridchan =               fread(filename, 1,  'int32'); % ??
header.specgridchannelselectval =   fread(filename, 1,  'int32');
header.specgriddatasize64 =         fread(filename, 1,  'int64'); % same as specgriddatasize
header.xstart =                     fread(filename, 1,  'int32');
header.xend =                       fread(filename, 1,  'int32');
header.ystart =                     fread(filename, 1,  'int32');
header.yend =                       fread(filename, 1,  'int32');
header.dummy =                      fread(filename, 234,'int32'); % skip the useless part of the header? Inconsistent with documentation
VZ =                                transpose(fread(filename, [3, header.vertpoints], 'single')); 
% VZ =                                fread(filename, [header.specgridchannelselectval, header.vertpoints], 'single'); % Broken: biases in every third entry of the first row
V =                                 single(VZ(:,1))/1000; % Broken
Z =                                 single(VZ(:,2)); % Broken
fullData =                           cell(header.yend, header.xend);

for y = header.ystart:header.yend
    for x = header.xstart:header.xend
        fullData{y, x} =             transpose(single(fread(filename, [header.specgridchan, header.vertpoints], 'single')));     
    end
end

fclose(filename);
% format(formatOrig);

%% Normalisations and massaging

% What I persume has something to do with the inverted Y-axis
fullData=flip(fullData,1);
channelData = zeros(header.xend, header.yend, header.vertpoints);

% Choose the channel of interest
for i=1:header.xend
    for j=1:header.yend
        channelData(i,j,:)=fullData{i,j}(:,options.channel);
    end
end

% Normalise
switch options.normalisation
    case 'raw'
        % Do nothing

    case 'median'
        for i=1:header.xend
            for j=1:header.yend
                channelData(i,j,:)=channelData(i,j,:)/(median(channelData(i,j,:)));   %median normalization          
        %         DATA(i,j,:)=DATA(i,j,:)/mean(abs(Data{i,j}(:,1)));                  %dIdV/average(abs(IV))
        %         IV(:)=Data{i,j}(:,1); temp(:)=DATA(i,j,:); DATA(i,j,:)=temp./IV;    %dIdV/IV
            end
        end
        
    case 'Imean'
        for i=1:header.xend
            for j=1:header.yend
%                 DATA(i,j,:)=DATA(i,j,:)/(median(DATA(i,j,:)));                      %median normalization          
                channelData(i,j,:)=channelData(i,j,:)/mean(abs(fullData{i,j}(:,1)));                  %dIdV/average(abs(IV))
        %         IV(:)=Data{i,j}(:,1); temp(:)=DATA(i,j,:); DATA(i,j,:)=temp./IV;    %dIdV/IV
            end
        end

    case 'IV'
        for i=1:header.xend
            for j=1:header.yend
%                 DATA(i,j,:)=DATA(i,j,:)/(median(DATA(i,j,:)));                      %median normalization          
        %         DATA(i,j,:)=DATA(i,j,:)/mean(abs(Data{i,j}(:,1)));                  %dIdV/average(abs(IV))
                IV(:)=fullData{i,j}(:,1); 
                temp(:)=channelData(i,j,:); 
                channelData(i,j,:)=temp./IV;    %dIdV/IV
            end
        end

    otherwise
        disp('Unknown normalisation setting: write your own or use the existing options!')
end

% Return to original numeric format for safety reasons
format(formatOrig.NumericFormat);
format(formatOrig.NumericFormat);


end