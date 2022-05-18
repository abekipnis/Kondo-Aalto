function [LineSpec] = CompLineFit(params, a, b, msOrig, E0Orig, MALine, ERange)
%LINEFITHANDLE Computes eigenmodes and produces a Line Spectra
%   Detailed explanation goes here

%% Initialisations
mCoeff = params(1);
global ms;
ms = mCoeff*ms;

E0Coeff = params(2);
global E0;
E0 = E0Coeff*E0;

EBroad = params(3).*1e-3;


%% Compute eigenmodes
[res, ~] = ComputeEigenmodes(a, b);

%% Compute and format line spectrum
NP = length(MALine);

LineSpec = ComputeLineSpectra(a, b, NP, res,...
    'EBroad', EBroad, ...
    'MALine', MALine,...
    'ERange', [min(ERange), max(ERange)],...
    'NE', length(ERange),...
    'PlotLine', false);

LineSpec = fliplr(LineSpec)./max(max(LineSpec));

%% Reset global parameters
ms = msOrig;
E0 = E0Orig;

end

