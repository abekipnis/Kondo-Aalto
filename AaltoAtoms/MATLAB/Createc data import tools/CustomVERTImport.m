function [Data, header] = CustomVERTImport(options)
    arguments
        options.files = 'no';
        options.path = 'no';
        options.normalisation string = 'IMax';
        options.direction string = 'all';
        options.BiasUnit char = 'mV'
        options.correctBiasOffset logical = true;
    end
    
%CUSTOMVERTIMPORT Imports all your Createc VERT files, however you like them
%   Imports VERT Files (see import functions for details), applies
%   normalisations.
%   Inputs (name-value pairs):
%       - files: the list of VERT files to be plotted. Default opens a file
%      browser.
%       - path: path to the folder containing the VERT files. Default opens
%       a file browser.
%       - normalisation: How to normalise the data
%           - 'raw' (default): plotting raw data (with NaN values removed)
%           - 'IMax': normalise all spectra with the current at the start
%           of the spectrum.
%           - 'StartNorm': normalise all spectra to the starting value
%           - 'avg': normalise each spectra to it's mean
%       - direction: which direction to take the spectra
%           - 'all' (default): plot all points in the files
%           - 'fwd': plot only the forward direction of bias sweeps
%           - 'bwd': -||- backward direction
%           - 'avg': averages the points in the forward and backward
%           directions
%       - BiasUnit: either 'mV' (default) or 'V' for volts
%       - correctBiasOffset: apply a simple bias offset (default true)
%   More optional inputs can be added as needed!
%   Outputs:
%       - h: handle to the generated figure
%       - data: A struct with a bunch of stuff like Bias, I, dIdV etc.
%       - header: Header of the first file to see what's what
%
% Adapted from PlotVERTFunction 3.9.2021
% Example on how to call the function:
%   >> % [Data, header] = CustomVERTImport('files', files, ...
%           'path', path,...
%           'direction', 'fwd',...
%           'normalisation', 'avg');

%% Import
% Use the supplied path and files if able
if strcmp(options.files, 'no')
    [files, path] = uigetfile('*.VERT',...
        'Select the files to plot', ...
        'Multiselect', 'on');
else
    files = options.files;
    path = options.path;
end

% Exit in case no files were imported
if isfloat(files)
    disp('No file selected: exiting')
    return;
end

%% Initialize necessary parameters
% Number of point spectra (NP). Also, read a VERT file to know what's what
if ischar(files)
    % One spectrum
    NP = 1;
    
    % Get the correct data line
    A = regexp(fileread([path files]),'\n','split');
    dataLine = find(contains(A,'DATA'), 1, 'last') + 2;
    
    [Bias, ~, ~] = ImportVERTfile([path files], [dataLine, inf]);
    header = ImportVERTHeader([path files]);

    if options.correctBiasOffset
        Bias = Bias - BiasOffset('file', files, 'path', path);
    end
else
    % Several spectra
    NP = length(files);
    
    % Get the correct data line
    A = regexp(fileread([path files{1}]),'\n','split');
    dataLine = find(contains(A,'DATA'), 1, 'last') + 2;
    
    [Bias, ~, ~] = ImportVERTfile([path files{1}], [dataLine, inf]);
    header = ImportVERTHeader([path files{1}]);

    % Line length from headers
    Pos1 = [header.XPos_nm, header.YPos_nm];

    Header2 = ImportVERTHeader([path files{end}]);
    Pos2 = [Header2.XPos_nm, Header2.YPos_nm];

    LineLength = norm(Pos2- Pos1);
%         LineDist = linspace(0, LineLength, NP);
    Line_nm = [0 LineLength];

    if options.correctBiasOffset
        Bias = Bias - BiasOffset('file', files{1}, 'path', path);
    end
end

% Get rid of NaN values
Bias = rmmissing(Bias);

% Remove duplicate points from the beginning and end
NSearch = 40;

[~,iStart,~] = unique(Bias(1:NSearch), 'last');
NEStart = iStart(end);

[~,iEnd,~] = unique(Bias(end-NSearch+1:end), 'first');
NEEnd = length(Bias) - NSearch + max(iEnd);

% Indexing works without sorting
BiasFilt = Bias(NEStart:NEEnd);

% Number of energy points: one way, both ways, or average
if strcmp(options.direction, 'all')
    NE = length(BiasFilt);
else
    NE = length(unique(BiasFilt));
end

% Mid point initialisation: reverts to NEEnd if the spectrum only has one
% direction
NEMid = find(Bias == min(Bias), 1);

% In case the starting bias is smaller than the end
if NEMid == NEStart
    NEMid = find(Bias == max(Bias) ,1);
end

% dIdV and I matrix initialisation
dIdV = zeros(NE, NP);
I = zeros(NE, NP);

%% Extract the data
for ii = 1:NP
    % Filename initialisation based on NP
    if NP ~= 1
        filename = [path files{ii}];
        VertGain = GetHeaderParam('VertManGain', 'files',files{ii}, 'path', path);
    else
        filename = [path files];
        VertGain = header.VertManGain;
    end
    
    [~, I_cur, LIX] = ImportVERTfile(filename, [dataLine, inf]);

    % Get rid of NaNs
    LIX = rmmissing(LIX);
    I_currm = rmmissing(I_cur);
    
    % Fix the current with preamp gain value
    % (Note that 'GainPreamplifier' is the scanning preamp gain)
    if VertGain ~= 9
        I_currm = I_currm.*10^(9 - VertGain);
    end

    %% Directions
    % If the spectra are taken in one direction, things get simple
    
    if length(BiasFilt) == length(unique(BiasFilt))
%         disp('No backward direction: plotting forward direction instead')
        dIdV(:,ii) = LIX(NEStart:NEEnd);
        I(:,ii) = I_currm(NEStart:NEEnd);

    else
         % Take the backward and forward channels
        LIXfw = LIX(NEStart:NEMid);
        LIXbw = flipud(LIX(NEMid:NEEnd));

        Ifw = I_currm(NEStart:NEMid);
        Ibw = flipud(I_currm(NEMid:NEEnd));
        
        switch options.direction
            case 'all'
                % won't work nicely with heat maps with several directions
                dIdV(:,ii) = LIX(NEStart:NEEnd);
                I(:,ii) = I_currm(NEStart:NEEnd);
                
            case 'avg'
                dIdV(:,ii) = (LIXfw + LIXbw)/2;.../mean(I_cur(1:NEStart));
                I(:,ii) = (Ifw + Ibw)/2;

            case 'fwd'
                dIdV(:,ii) = LIXfw;
                I(:,ii) = Ifw;

            case 'bwd'
                % If the spectra are only taken in one direction, take the
                % forward direction instead

                dIdV(:,ii) = LIXbw;
                I(:,ii) = Ibw;

        end
    end
    
    %% Normalisations
    switch options.normalisation
        case 'raw'
            continue
            
        case 'IMax'
            dIdV(:,ii) = dIdV(:,ii)./I(1, ii);
            
        case 'StartNorm'
            dIdV(:,ii) = dIdV(:,ii)./dIdV(1, ii);
            
        case 'avg'
            dIdV(:,ii) = dIdV(:,ii)./mean(dIdV(:,ii));
    end

end

%% Adjust the bias vector based on directions
switch options.direction
    case 'all'
        % Default initialisation is good: no need to do anything
    case 'bwd'
        BiasFilt = flipud(Bias(NEMid:NEEnd));
    otherwise
        BiasFilt = Bias(NEStart:NEMid);
end

%% Bias unit
if strcmp(options.BiasUnit, 'V')
    BiasFilt = 1e-3.*BiasFilt;
end

%% Produce the data struct
% Always available ones
Data.files = files;
Data.path = path;
Data.dIdV = dIdV;
Data.Bias = BiasFilt;
Data.I = I;
Data.NP = NP;
Data.NE = NE;

% Conditional ones
if NP > 1
    Data.LineLength = LineLength;
else
    Data.LineLength = [];
end


end

