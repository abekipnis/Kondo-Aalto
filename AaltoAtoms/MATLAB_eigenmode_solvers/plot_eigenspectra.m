InitializeGlobals('Ag');

radii = 2:0.01:5;
eigs = [];
for r = radii
    disp(r);
    [res_vac, model_vac] = ComputeEigenmodes(r, r, ...
    'plotAll', 0,...
    'HMax', 0.1,...
    'energyRange',300e-3);
    eigs = [eigs res_vac.Eigenvalues(1)];
end

T = table(radii, eigs+E0, 'VariableNames',{'radii','eigenvals (eV)'});
writetable(T,'Occupied_eigvals.txt')



