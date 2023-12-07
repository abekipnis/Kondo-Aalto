function [param] = GetHeaderParam(paramName, options)
    arguments
        paramName char
        options.files = 'no'
        options.path = 'no'
    end

%GETHEADERPARAM Returns the value(s) of a VERT file header parameter
% Inputs:
%   - paramName: char with the parameter name needed (see ImportVertHeader
%   for the list of parameter names)
%   - files: the files from which the header will be extracted. Default
%   opens a file browser.
%   - path: path of said files.
% Output:
%   - param: a vector of numerical parameter values
%
% FIXME: One path per function call: use loops to go through several paths
% FIXME: Make param a struct or cell array to store multiple types of
% header parameters.

%% File names and path
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

%% Do the header imports in a loop
% Number of point spectra (NP)
if ischar(files)
    % One spectrum
    NP = 1;
else
    % Several spectra
    NP = length(files);
end

% FIXME: taking only numerical values for now
param = zeros(1,NP);

if NP == 1
    header = ImportVERTHeader([path '/' files]);
%     param = getfield(header, paramName);
    param = header.(paramName);
else
    for ii = 1:NP
        header = ImportVERTHeader([path '/' files{ii}]);
        param(ii) = header.(paramName);
    end
end

