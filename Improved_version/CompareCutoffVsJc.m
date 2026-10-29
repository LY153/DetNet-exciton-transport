% CompareCutoffVsJc.m
% Sweeps the cutoff parameter Jc and plots, versus x:
%   left y-axis  - mean(|dE|/|E0|), the relative eigenvalue shift
%   right y-axis - mean(1 - overlap_corrected), the eigenstate deviation
% The figure is saved as CutoffCompareResult.fig next to this script.

clear; clc; close all;

m       = 5;
xlist   = 4:8;
Jc_list = [1 5 10];

nX = numel(xlist);
nJ = numel(Jc_list);

mean_rel_dE = zeros(nX, nJ);
mean_dev    = zeros(nX, nJ);

for i = 1:nX
    for j = 1:nJ
        r = CompareHamiltonians(m, xlist(i), 'Jc', Jc_list(j), 'verbose', false);
        mean_rel_dE(i,j) = mean(abs(r.dE) ./ abs(r.E0));
        mean_dev(i,j)    = mean(r.dev_corrected);
        fprintf('x=%d, Jc=%-4g : mean|dE|/E0 = %.4e, mean(1-overlap_corrected) = %.4e\n', ...
            xlist(i), Jc_list(j), mean_rel_dE(i,j), mean_dev(i,j));
    end
end

fig = figure('Name', 'CutoffCompareResult', 'Color', 'w', 'Position', [100 100 760 520]);

colors = get(groot, 'defaultAxesColorOrder');
colors = colors(1:nJ, :);

yyaxis left
hL = gobjects(1, nJ);

hL(1) = plot(xlist, mean_rel_dE(:,1), 'o-', 'LineWidth', 1.5, ...
    'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
hold on;
hL(2) = plot(xlist, mean_rel_dE(:,2), '^-', 'LineWidth', 1.5, ...
    'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
hold on;
hL(3) = plot(xlist, mean_rel_dE(:,3), 's-', 'LineWidth', 1.5, ...
    'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
hold on;

ylabel('$\overline{|\Delta E|/E_0}$','Interpreter','latex');
hold off;

yyaxis right
hR = gobjects(1, nJ);

hR(1) = plot(xlist, mean_dev(:,1), 'o--', 'LineWidth', 1.5, ...
    'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
hold on;
hR(2) = plot(xlist, mean_dev(:,2), '^--', 'LineWidth', 1.5, ...
    'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
hold on;
hR(3) = plot(xlist, mean_dev(:,3), 's--', 'LineWidth', 1.5, ...
    'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
hold on;

ylabel('$\overline{(1-|\langle\psi_0|\psi\rangle|^2)}$','Interpreter','latex');
hold off;

xlabel('$x$','Interpreter','latex');
grid on;
xticks(xlist);

legLabels = cell(1, 2*nJ);
for j = 1:nJ
    legLabels{j}    = ['$J_c=' num2str(Jc_list(j)) '$ : $\overline{|\Delta E|/E_0}$'];
    legLabels{nJ+j} = ['$J_c=' num2str(Jc_list(j)) '$ : $\overline{(1-|\langle\psi_0|\psi\rangle|^2)}$'];
end
legend([hL, hR], legLabels, 'Location', 'northwest', 'Interpreter', 'latex');

out_file = fullfile(fileparts(mfilename('fullpath')), 'CutoffCompareResult.fig');
savefig(fig, out_file);
fprintf('Saved %s\n', out_file);
