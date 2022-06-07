clear all
close all
clc

% Plots VERT files: plots differently based on the number of files given
%   - If only one file is supplied, the spectrum is plotted normally
%   - If several spectra are opened, they are plotted as a heat map
%   (interpreted as line spectra)

%% Read the files from a folder
% Vert to heat map
% Dum version: files in the same directory

% list files
% list = ls('*.VERT');

% Less dum version: choose the files to open
[files, path] = uigetfile('*.VERT',...
    'Select the files to plot', ...
    'Multiselect', 'on');
% 
% For easy access, save the path and file to a variable and load it here
% files = 'Createc2_201223.191235.L0050.VERT';
% path = '\\home.org.aalto.fi\aaprom1\data\Documents\Väikkäri\Kokeet\Coupling of magnetic atoms in quantum corrals\Results and analysis\Relevant data\Corral C (23.12.)\';


%% Case: only one vert file opened
NP = length(files);

if isfloat(files)
    disp('No file selected: exiting')
else

    if ischar(files)
        [Bias, I, LIX] = ImportVERTfile([path files]);

        % Plot
        figure;
        plot(Bias, LIX, '-o')
        title(files, ...
            'Interpreter', 'none')
        xlabel('Bias, mV')
        ylabel('dI/dV, a.u.')
        set(gca, 'FontSize', 14)
    else
        %% Case 2: several files chosen (a line spec)
        % Read one vert file to know what's what
        [Bias, ~, ~] = ImportVERTfile([path files{1}]);
        
        % Get rid of NaN values
        Bias = rmmissing(Bias);
        Bias = Bias(2:end);
        NE = round(length(Bias));
        % NP = size(list, 1);
        NP = length(files);
        dIdV = zeros(NE, NP);
        I = zeros(NE, NP);
        
        % Line length from headers
        Header1 = ImportVERTHeader([path files{1}]);
        Pos1 = [Header1.XPos_nm, Header1.YPos_nm];

        Header2 = ImportVERTHeader([path files{end}]);
        Pos2 = [Header2.XPos_nm, Header2.YPos_nm];

        LineLength = norm(Pos2- Pos1);
        LineDist = linspace(0, LineLength, NP);
        Line_nm = [0 LineLength];
        %% If only one Vert file is chosen, plot straight away
        % if NP == 1
        %     

        %% for list 
        % matrix(ii,res) = ImportVERTfile(list(ii))

        for ii = 1:NP
            filename = [path files{ii}];
            [~, I_cur, LIX] = ImportVERTfile(filename);
            
            % Get rid of NaNs
            LIX = rmmissing(LIX);
            I_cur = rmmissing(I_cur);
            
            dIdV(:,ii) = LIX(1:NE);%./I(1);
            I(:,ii) = I_cur(1:NE);
        end

        %% Plot heat map
        figure
        imagesc(Line_nm, [Bias(1) Bias(NE)], dIdV)
        set(gca, 'YDir', 'normal')
        xlabel('Distance, nm')
        ylabel('Bias, mV')
        set(gca, 'FontSize', 14)
        c = colorbar;
        set(c, 'YTickLabel', [])
        ylabel(c, {'{\it dI/dV}, a.u.'})
        
    end
end

%% Plot Lorentzian
% 
% Broadening = @(E, Ec, width) 1./(pi.*width).*...
%     (width.^2./((E - Ec).^2 + width.^2));
% 
% E_cent = -30.49;
% width = 6.5;
% NE = round(length(Bias)/2);
% amp = 2300;
% 
% hold on
% % figure;
% plot(Bias, amp.*Broadening(Bias, E_cent, width) + median(LIX(3:NE)), '-r')