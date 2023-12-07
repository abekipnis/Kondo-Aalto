function [offset] = BiasOffset(options)

arguments
    options.file char = 'no'
    options.path char = 'no'
end

%BIASOFFSET Check the bias offset of a VERT spectrum
%   Call in a loop for several spectra
%   Optional name-value pair inputs
%       - file: VERT file name
%       - path: path to said VERT file
%   Output:
%       - offset: the bias where the current crosses zero [mV]
%
%   NOTE: USE WITH METALS ONLY! Does not work with a gap at Fermi

%% Import the data
% Get a file and path in case it wasn't supplied
if strcmp(options.file, 'no')
    [file, path] = uigetfile('*.VERT',...
        'Select the files to plot', ...
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

% Get the correct data line
A = regexp(fileread([path '\' file]),'\n','split');
dataLine = find(contains(A,'DATA'), 1, 'last') + 2;

[Bias, I, ~] = ImportVERTfile([path '\' file], [dataLine, inf]);

%% Find the bias where current crosses zero

% If the spectra have several directions, use only the forward direction
[Bias, ia, ~] = unique(Bias, 'stable');
I = I(ia);

% Also, make sure the current vector is unique. Sometimes it is not.
[I, ib, ~] = unique(I, 'stable');
Bias = Bias(ib);

% Interpolate to get the offset
offset = interp1(I, Bias, 0);


end

