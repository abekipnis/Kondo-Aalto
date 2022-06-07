clear all
close all
clc

%% TEST FUNCTION SOLVERS

InitializeGlobals('Ag');

global ms

ms = ms/1.1;

a_vac = 3.5;
b_vac = 2.65;

[res_vac, model_vac] = ComputeEigenmodes(a_vac, b_vac, ...
    'plotAll', 0,...
    'HMax', 0.1);

NP = 100;

LSp = ComputeLineSpectra(a_vac, b_vac, NP, res_vac,...
    ...'YLim', [-15 45],...
    'PlotLine', true);

ComputePointSpectrum(a_vac, b_vac, res_vac, model_vac, [0,0])