function [Bias, I, LIX] = ImportVERTfile(filename, dataLines)
%IMPORTVERTFILE Import data from a VERT file
%  [Bias, I, LIX] =
%  IMPORTFILE(FILENAME) reads data from VERT file FILENAME for the
%  default selection.  Returns the data as column vectors.
%
%  [Bias, I, LIX] =
%  IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  [Bias, I, LIX] = ImportVERTfile("C:\LocalUserData\User-data\aaprom1\Matlab\Import functions\Createc2200306.153317.L0008.VERT", [617, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 20-Oct-2020 13:23:36


%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [630, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ["\t", " ", "="];

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3"], "ThousandsSeparator", ",");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
Bias = tbl.VarName2;

% Current units: search DAC from STMAFM help file for conversion 
% HARDCODED GAIN=9! Switch in the calling function/script based on header
I = tbl.VarName5.*10/2^19/10^9;
LIX = tbl.VarName6;


end