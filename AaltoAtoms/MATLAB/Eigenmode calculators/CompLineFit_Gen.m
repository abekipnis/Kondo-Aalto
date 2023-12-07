function [LineSpec] = CompLineFit_Gen(options)
arguments
    options.a double
    options.b double
    options.msOrig double = 2.27424e-12;    % Default: Ag(111)
    options.E0Orig double = -0.067;         % Default: Ag(111)
    options.ERange
    options.MALine
    options.meffCoeff double = 1;
    options.E0Coeff double = 1;
    options.EBroad double = 5;
    options.CenterAtom logical = false;
    options.FocusAtom logical = false;
    options.centerPc logical = false;
    options.atomPotential double = 0.9;
    options.atomRadius double = 0.6;
    options.WeissDistance double = 3;
    options.DecayFactor double = 0.5;
    options.StepHeight double = 0;
    options.plotAll logical = false;
    options.HMax double = 0.3;
    options.plotLine logical = false;
    options.normaliseLine logical = true;
    options.MALineOffset double = 0;
    options.EExtension double = 0;
end

%COMPLINEFIT_GEN Computes eigenmodes and produces line spectra
%   Make an anonymous handle to this function to produce fits with specific
%   parameters
%   Inputs:
%       - a: major axis
%       - b: minor axis
%       - msOrig: original effective mass value (to bypass problems with
%       global parameters)
%       - E0Orig: -||- for surface state onset
%       - ERange: Energy range vector to compute the line spectra, [eV]
%       - MALine: points along the major axis, corral center at 0 [nm]
%       ...
%       (see the parameter descriptions from ComputeEigenmodes and
%       ComputeLinespectra)
%       ...
%       - normaliseLine: flag for normalising the line spectrum to it's
%       maximum. Set this to true while fitting!
%       - MALineOffset: Moves the major axis line by a constant, [nm]
%       - EExtension: Adds a constant amount of energy for eigenmode
%       calculations wrt. max(ERange), [eV]
% FIXME:
%   - Default WeissDistance will cause issues if it isn't calculated prior
%   to function call

%% Initialisations
% Global parameters removed to make the functions run in parallel for-loops
% mCoeff = options.meffCoeff;
% global ms;
% ms = mCoeff*ms;
% 
% E0Coeff = options.E0Coeff;
% global E0;
% E0 = E0Coeff*E0;

EBroad = options.EBroad.*1e-3;


%% Compute eigenmodes
[res, ~] = ComputeEigenmodes(options.a, options.b,...
    'atomPotential', options.atomPotential,...
    'energyRange', max(options.ERange) + options.EExtension,...
    'HMax',options.HMax,...
    'atomRadius',options.atomRadius,...
    'centerAtom',options.CenterAtom, ...
    'focusAtom', options.FocusAtom,...
    'centerPc',options.centerPc,...
    'plotAll', options.plotAll, ...
    'ms_local', options.msOrig*options.meffCoeff,...
    'E0_local', options.E0Orig*options.E0Coeff);

%% Compute and format line spectrum
NP = length(options.MALine);

LineSpec = ComputeLineSpectra(options.a, options.b, NP, res,...
    'EBroad', EBroad, ...
    'MALine', options.MALine - options.MALineOffset,...
    'ERange', [min(options.ERange), max(options.ERange)],...
    'NE', length(options.ERange),...
    'PlotLine', options.plotLine,...
    'DecayFactor',options.DecayFactor,...
    'WeissDistance',options.WeissDistance,....
    'ms_local', options.msOrig*options.meffCoeff,...
    'E0_local', options.E0Orig*options.E0Coeff);

% Normalise to maximum value: can make it seem that different energy ranges
% produce different eigenmode and point spectrum results. 
if options.normaliseLine
    LineSpec = fliplr(LineSpec)./max(max(LineSpec));
end

% Add a surface state onset step
FD = @(A, E, Ec, w) A - 1.*A.*(exp((E-Ec)./w) + 1).^-1;
% step = FD(options.StepHeight, options.ERange, E0, 5e-3);
step = FD(options.StepHeight, options.ERange, options.E0Orig*options.E0Coeff, 5e-3);
LineSpec = LineSpec + step';

% %% Reset global parameters
% ms = options.msOrig;
% E0 = options.E0Orig;

end

