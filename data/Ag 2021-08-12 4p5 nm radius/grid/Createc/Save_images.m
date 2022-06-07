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

DATA = zeros(header.xend, header.yend, header.vertpoints);

for i=1:header.xend
    for j=1:header.yend
        if (sum(size(Data{i,j})==0)==2 || size(Data{i,j}(:,5),1)~=header.vertpoints)
            DATA(i,j,:)=zeros(1,header.vertpoints);
        else 
            DATA(i,j,:)=Data{i,j}(:,6);
        end
    end
end

a=size(DATA);
blx=linspace(0,header.dx,a(1));  
bly=linspace(0,header.dy,a(2)); 

%__________________________________________________________________________
%ANALYSIS

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

figure(1)
set(gcf, 'Position',  [100, 100, 500, 500])
for i=1:a(3)
    cond_map(1:a(1),1:a(2))=DATA(:,:,i);
    surf(blx, bly, cond_map, 'Edgecolor', 'none');
        
    vec=sort(reshape(cond_map,1,a(1)*a(2)));        
    caxis([vec(round(a(1)*a(2)*0.05)) vec(round(a(1)*a(2)*0.95))])
    axis equal
    xlim([blx(1) blx(end)])
    ylim([bly(1) bly(end)])
    view([0 0 1])
    
    saveas(gcf,['dIdV_map_at_', num2str(xdata(i)),'V.png'])
end



