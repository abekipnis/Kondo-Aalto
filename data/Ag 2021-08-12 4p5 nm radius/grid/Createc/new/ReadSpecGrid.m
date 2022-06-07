function [header, V, Z, spectra] = ReadSpecGrid(file)

%   This function reads the specgrid binary file.
%   Inputs
%       filename: .specgrid file only
%   Outputs
%       header: struct
%       V,Z:    column matrix
%       data:   cell array

%   HEADER
%   version                     integer
%   nx,ny                       integer
%   dx,dy                       integer
%   specxgrid,specygrid         integer
%   vertpoints                  integer
%   vertmandelay                integer
%   vertmangain                 integer
%   biasvoltage                 single
%   tunnelcurrent               single
%   imagedatasize               integer
%   specgriddatasize            integer
%   specgridchan                integer
%   specgridchannelselectval    integer
%   specgriddatasize64          int64
%   xstart,xend,ystart,yend     integer
%   specgridchannelselectval2   integer
%   specgridnx  		integer
%   specgridny 			integer
%   specgriddx 			integer
%   specgriddy 			integer
%   specgridcenterx 		integer
%   specgridcentery             integer
%   dummy                       array[1..227] of integer

%   V AND Z
%   V, Z                        single
%   Length = vertpoints 

%   SPECTRA
%   spectra                     single
%   Spectrum (X=1, Y=1), 
%       Length of spectrum per grid point: specgridchan*vertpoints 
%   Spectrum (X=2, Y=1), 
%   .....
%   Spectrum (X=xend, Y=yend)

%   Conversion table
%   Delphi      Matlab
%   integer     int32
%   single      single
%   int64       int64

filename =                   fopen(file);
if filename == -1
    error('File not found or can not be read!')
end

header.ver =                        fread(filename, 1,  'int32');
header.nx =                         fread(filename, 1,  'int32');
header.ny =                         fread(filename, 1,  'int32');
header.dx =                         fread(filename, 1,  'int32');
header.dy =                         fread(filename, 1,  'int32');
header.specxgrid =                  fread(filename, 1,  'int32');
header.specygrid =                  fread(filename, 1,  'int32');
header.vertpoints =                 fread(filename, 1,  'int32');
header.vertmandelay =               fread(filename, 1,  'int32');
header.vertmangain =                fread(filename, 1,  'int32');
header.biasvoltage =                fread(filename, 1,  'single');
header.tunnelcurrent =              fread(filename, 1,  'single');
header.imagedatasize =              fread(filename, 1,  'int32');
header.specgriddatasize =           fread(filename, 1,  'int32');
header.specgridchan =               fread(filename, 1,  'int32');
header.specgridchannelselectval =   fread(filename, 1,  'int32');
header.specgriddatasize64 =         fread(filename, 1,  'int64');
header.xstart =                     fread(filename, 1,  'int32');
header.xend =                       fread(filename, 1,  'int32');
header.ystart =                     fread(filename, 1,  'int32');
header.yend =                       fread(filename, 1,  'int32');
header.specgridchannelselectval2=   fread(filename, 1,  'int32');
header.specgridnx =                 fread(filename, 1,  'int32');
header.specgridny =                 fread(filename, 1,  'int32');
header.specgriddx =                 fread(filename, 1,  'int32');
header.specgriddy =                 fread(filename, 1,  'int32');
header.specgridcenterx =            fread(filename, 1,  'int32');
header.specgridcentery =            fread(filename, 1,  'int32');
header.dummy=                       fread(filename,227, 'int32'); 

VZ =                                transpose(fread(filename, [header.specgridchannelselectval, header.vertpoints], 'single'));
V =                                 single(VZ(:,1));
Z =                                 single(VZ(:,2));
spectra =                           cell(header.yend, header.xend);

for y = 1:header.specgridny
    for x = 1:header.specgridnx
        spectra{y, x} =             transpose(single(fread(filename, [header.specgridchan, header.vertpoints], 'single')));     
    end
end

fclose(filename);

end


