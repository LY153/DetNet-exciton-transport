% run_cutoff_compare.m
% Runs CompareHamiltonians for m = 5, x = 4..8 and saves the results to
% CutoffCompare.mat next to this script.

m = 5;
xlist = 4:8;

CutoffCompare = struct();
CutoffCompare.m = m;
CutoffCompare.x = xlist(:).';
CutoffCompare.combos = cell(numel(xlist), 1);

for k = 1:numel(xlist)
    r = CompareHamiltonians(m, xlist(k));
    if istable(r.corrected)
        r.corrected = table2array(r.corrected);
    end
    CutoffCompare.combos{k} = r;
end

out_file = fullfile(fileparts(mfilename('fullpath')), 'CutoffCompare.mat');
save(out_file, 'CutoffCompare');
fprintf('Saved %s\n', out_file);
