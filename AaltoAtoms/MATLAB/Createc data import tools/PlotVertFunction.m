function [h, Data, header] = PlotVertFunction(options)
    arguments
        options.files = 'no';
        options.path = 'no';
        options.correctBiasOffset logical = true;
        options.normalisation string = 'IMax';
        options.BGSpectrumFile = 'no';
        options.BGSpectrumPath = 'no';
        options.direction string = 'all';
        options.multiPlotType = 'heatMap';
        options.stackDisp = 0.1;
        options.stackColorMap char = 'cool';
        options.xAxisType = 'dist'
        options.customXLimits double = [0, 1];
        options.customXLabel = 'Custom Label'
        options.BiasUnit char = 'mV'
    end
    
%PLOTVERTFUNCTION Plots all your Createc VERT files
%   Imports VERT Files (see import functions for details), applies
%   normalisations, and plots the spectra in multiple ways.
%   Inputs (all optional name-value pairs):
%       - files: the list of VERT files to be plotted (char vector/string 
%       for individual spectrum, cell array for multiple spectra). Default 
%       'no' opens a file browser.
%       - path: path to the folder containing the VERT files (char vector). 
%       Default 'no' opens a file browser.
%       - correctBiasOffset (default true): apply a simple bias offset (see
%       function BiasOffset for details)
%       - normalisation: How to normalise the data
%           - 'raw' (default): plotting raw data (with NaN values removed)
%           - 'IMax': normalise all spectra with the current at the start
%           of the spectrum.
%           - 'StartNorm': normalise all spectra to the starting LIX value
%           - 'avg': normalise each spectrum to its mean value
%       - BGSpectrumFile: name of a VERT file subtracted from all the data,
%       after applying the specified normalisation
%       - BGSpectrumPath: path to the background VERT spectrum
%       - direction: which direction to take the spectra
%           - 'all' (default): plot all points in the files
%           - 'fwd': plot only the forward direction of bias sweeps
%           - 'bwd': -||- backward direction
%           - 'avg': averages the points in the forward and backward
%           directions
%       - multiPlotType: how to plot multiple spectra
%           - 'heatMap' (default): plot the spectra as a heat map
%           - 'stack': plots the curves individually with a vertical
%           shift between curves
%       - stackDisp: a constant displacement of the curves for a stacked
%       curve plot (remember to adjust according to normalisation)
%       - stackColorMap: a char vector specifying a color map for stack
%       plots (default 'cool')
%       - xAxisType: What kind of x-axis to use in heatmaps
%           - 'n': Spectrum index, from the sequence in the file list
%           - 'dist' (default): position coordinate differences, for line spectra
%           - 'custom': User-defined variable and label.
%       - customXLimits: vector [xmin, xmax] where xmin < xmax. Only works
%       if xAxisType is set to 'custom'!
%       - customXLabel: char vector with custom x label. Only works if 
%       xAxisType is set to 'custom'!
%       - BiasUnit: either 'mV' (default) or 'V' for volts
%   More optional inputs can be added as needed!
%
%   Outputs:
%       - h: handle to the generated figure
%       - data: A struct with data and parameters extracted from spectra
%       - header: Header of the first spectrum file
%
% Example on how to call the function:
%   >>  [h1, Data, header] = PlotVertFunction('files', files, ...
%           'path', path,...
%           'direction', 'fwd',...
%           'multiPlotType', 'stack',...
%           'normalisation', 'avg',...
%           'stackDisp', 0.03);
% TO DO:
% - Optional saving of the file path
% - Specifying customXLabel or customXLimits should generate custom axes 
%   regardless of xAxisType
% - All imported files use the same bias axis: reasonable for line spectra 
%   but not so much for other uses

%% Import
% Use the supplied path and files if able
if strcmp(options.files, 'no') || strcmp(options.path, 'no')
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
% Number of point spectra (NP). Read the first VERT file to know what's what
if ischar(files)
    % One spectrum
    NP = 1;
    
    % Get the correct data line
    A = regexp(fileread([path files]),'\n','split');
    dataLine = find(contains(A,'DATA'), 1, 'last') + 2;
    
    % Import the bias voltage vector
    [Bias, ~, ~] = ImportVERTfile([path files], [dataLine, inf]);
    header = ImportVERTHeader([path files]);
    
    % Correct for bias offset
    if options.correctBiasOffset
        Bias = Bias - BiasOffset('file', files, 'path', path);
    end
else
    % Several spectra
    NP = length(files);
    
    % Get the correct data line
    A = regexp(fileread([path files{1}]),'\n','split');
    dataLine = find(contains(A,'DATA'), 1, 'last') + 2;

    % Import the bias voltage vector
    [Bias, ~, ~] = ImportVERTfile([path files{1}], [dataLine, inf]);
    header = ImportVERTHeader([path files{1}]);

    % Line length from headers
    Pos1 = [header.XPos_nm, header.YPos_nm];
    
    % Import the header of the last file in the list
    Header2 = ImportVERTHeader([path files{end}]);
    Pos2 = [Header2.XPos_nm, Header2.YPos_nm];

    LineLength = norm(Pos2- Pos1);
    Line_nm = [0 LineLength];

    % Correct for bias offset
    if options.correctBiasOffset
        % FIXME: make fancier (mean or point-wise)
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
    NEMid = find(Bias == max(Bias), 1);
end

% dIdV and I matrix initialisation
dIdV = zeros(NE, NP);
I = zeros(NE, NP);

%% Extract the data

for ii = 1:NP
    % Filename and vertGain initialisation based on NP
    if NP ~= 1
        filename = [path files{ii}];
        VertGain = GetHeaderParam('VertManGain', 'files',files{ii}, 'path', path);
    else
        filename = [path files];
        VertGain = header.VertManGain;
    end
    
    [~, I_cur, LIX] = ImportVERTfile(filename, [dataLine, inf]);

    % Get rid of NaN values, if any
    LIX = rmmissing(LIX);
    I_currm = rmmissing(I_cur);

    % Fix the current (and dIdV) with preamp gain value
    if VertGain ~= 9
        I_currm = I_currm.*10^(9 - VertGain);
        LIX = LIX.*10^(9 - VertGain);
    end

    %% Directions
    % If the spectra are taken in one direction, things get simple
    
    % FIXME: BiasFilt doesn't update in the loop: non-uniform data gets
    % messed up regardless of the plot type
    if length(BiasFilt) == length(unique(BiasFilt))
%         disp('No backward direction: plotting forward direction instead')
%         disp(num2str(ii));
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
    % Normalise the dI/dV data
    switch options.normalisation
        case 'raw'
            continue
            
        case 'IMax'
            % First current value
            dIdV(:,ii) = dIdV(:,ii)./I(1, ii);
            
        case 'StartNorm'
            % First dIdV value
            dIdV(:,ii) = dIdV(:,ii)./dIdV(1, ii);
            
        case 'avg'
            % Average dI/dV value
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

%% Normalize by a background spectrum
% FIXME: BG spectrum must have identical number of points and bias range
% FIXME: Only one file can be selected
if ~strcmp('no', options.BGSpectrumFile)
    if strcmp('no', options.BGSpectrumPath)
        options.BGSpectrumPath = path;
    end

    [BGData, ~] = CustomVERTImport("files", options.BGSpectrumFile, ...
        "path", options.BGSpectrumPath, ...
        "normalisation", options.normalisation);

    % Subtraction (division only works for properly normalised data)
    dIdV = dIdV - BGData.dIdV;
end

%% Plot the spectra

h = figure;

% Adjust bias unit: default unit in Createc output files is mV
if strcmp(options.BiasUnit, 'V')
    BiasFilt = 1e-3.*BiasFilt;
    BiasLabel = 'Bias (V)';
else
    BiasLabel = 'Bias (mV)';
end

% Start plotting, type based on NP
if NP ==1
    % One spectrum: classic line plot
    plot(BiasFilt, dIdV, '-o')
    title(files, ...
        'Interpreter', 'none')
    xlabel(BiasLabel)
    ylabel('d{\itI}/d{\itV} (arb. units)')
    set(gca, 'FontSize', 14)
else
    switch options.multiPlotType

        case 'heatMap'
            % Set the x-axis
            switch options.xAxisType
                case 'dist'
                    xAxis = Line_nm;
                    xLabel = 'Distance (nm)';
                    
                case 'n'
                    xAxis = [1 NP];
                    xLabel = 'Point';
                    
                case 'custom'
                    xAxis = options.customXLimits;
                    xLabel = options.customXLabel;
            end
            
            imagesc(xAxis, [BiasFilt(1) BiasFilt(end)], dIdV)
            set(gca, 'YDir', 'normal')
            xlabel(xLabel)
            ylabel(BiasLabel)
            set(gca, 'FontSize', 14)
            c = colorbar;
            set(c, 'YTickLabel', [])
            ylabel(c, {'d{\itI}/d{\itV} (arb. units)'})
            
        case 'stack'
            % Plot color setup
            fullColors = colormap(options.stackColorMap);
            
            hold on
            for jj = 1:NP
                plot(BiasFilt, dIdV(:,jj) + options.stackDisp*jj, ...
                    "Color", fullColors(round(jj*256/NP), :))
            end
            axis tight
            xlabel(BiasLabel)
            ylabel('d{\itI}/d{\itV} (arb. units)')

            % Set the x-axis
            switch options.xAxisType
                case 'dist'
                    TL = Line_nm;
                    cLabel = 'Distance (nm)';

                case 'n'
                    TL = [1 NP];
                    cLabel = 'Point';

                case 'custom'
                    TL = options.customXLimits;
                    cLabel = options.customXLabel;
            end

            % Mess with the colorbar: adjust decimal precision as needed
            c = colorbar('Ticks', [0 1], 'TickLabels', arrayfun(@(x) sprintf('%.3f',x),TL,'un',0));
%             c = colorbar('Ticks', linspace(0,1, NP), 'TickLabels', linspace(1,NP, NP));
            c.Label.String = cLabel;
            c.Label.VerticalAlignment = 'bottom';
            
            set(gca, 'FontSize', 14)
    end
end

%% Produce the Data struct

Data.files = files;
Data.path = path;
Data.dIdV = dIdV;
Data.Bias = BiasFilt;
Data.I = I;
Data.NP = NP;
Data.NE = NE;

% Optional content for multiple spectra
if NP == 1
    Data.LineLength = 0;
else
    Data.LineLength = LineLength;
end

end

