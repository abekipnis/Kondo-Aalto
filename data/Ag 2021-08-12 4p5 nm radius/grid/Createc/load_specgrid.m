clear all; close all; clc;
format longE

%__________________________________________________________________________
%LOADING DATA

disp( 'choose data directory' )
[baseName, folder] = uigetfile('*.specgrid');
[header, xdata, Z, Data] = ReadSpecGrid([folder baseName]); 
mkdir([folder 'Analysis\'])
disp( 'loading data...' )
Data=flip(Data,1);

for i=1:header.xend
    for j=1:header.yend
%          if (sum(size(Data{i,j})==0)==2)% || size(Data{i,j}(:,6),1)~=header.vertpoints)
%             DATA(i,j,:)=zeros(1,header.vertpoints);
%         else 
            DATA(i,j,:)=Data{i,j}(:,6);
%         end
    end
end

a=size(DATA);
blx=linspace(0,header.dx,a(1));    %creating [nm] vector based on window size and number of points
bly=linspace(0,header.dy,a(2));    %creating [nm] vector based on window size and number of points
% xdata=xdata/100;    

%__________________________________________________________________________
%ANALYSIS

dI_su=0;
for i=1:header.xend
    for j=1:header.yend
%         DATA(i,j,:)=DATA(i,j,:)/(median(DATA(i,j,:)));                      %median normalization          
%         DATA(i,j,:)=DATA(i,j,:)/mean(abs(Data{i,j}(:,1)));                  %dIdV/average(abs(IV))
         try
            IV(:)=Data{i,j}(:,1);
            temp(:)=DATA(i,j,:); 
            DATA(i,j,:)=temp./IV;    %dIdV/IV
        catch 
            DATA(i,j,:)= DATA(i,j,:);
        end
    end
end
dI_sum(:)=median(median(DATA));

QPI_movie(xdata, dI_sum, DATA, folder, baseName, blx, bly)