function [files, path] = SavePath(dataName)
%SAVEPATH Saves filenames and path into a .MAT file
%   Input:
%       - dataName: string used to name the files (without file extension)
%   Output: 
%       - files: cell array or string with the file name(s)
%       - path: path string terminated with a '\'

 [files, path] = uigetfile('*.VERT',...
        'Select the files to plot', ...
        'Multiselect', 'on');

 save([dataName '.mat'], 'path', 'files');

end

