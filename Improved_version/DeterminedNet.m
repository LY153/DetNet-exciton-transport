% DeterminedNet  Optimized DetNet driver: runs the simulation for m = m_f
% layers and x = x_r, then plots the site-averaged efficiency eta vs x and vs NH.
clear; close all;

m_f = 5;
x_r = 4:8;
T_max = 1e7;
convergence_threshold = 0.99;
save_interval = 500;
resume_mode   = false;             % true = resume computation, false = restart computation

numCombos = length(x_r);

data = zeros(numCombos, 4);
results_raw = cell(numCombos, 2);
results_avg = cell(numCombos, 2);
eta_details = cell(numCombos, 1);

delete(gcp('nocreate'));
try
    n_cores = feature('numcores');
catch
    n_cores = 4;
end
n_workers = max(1, min(12, n_cores));   % raise 12 if you have more cores
parpool(n_workers);

fprintf('m=%d, x=%s, T_max=%d\n', m_f, mat2str(x_r), T_max);

tic_all = tic;
for idx = 1:numCombos
    x = x_r(idx);
    m = m_f;

    first_site = m*(m-1)*x/2 + 2;
    last_site  = m*(m+1)*x/2 + 1;
    site_indices = first_site : last_site;
    M_sites = length(site_indices);

    state_file = fullfile(tempdir, sprintf('state_m%d_x%d.mat', m, x));

    [eta, t_elapsed, population_s_cell, totalpopulation_cell, ~] = ...
        DetNetcalculation(m, x, site_indices, ...
            T_max, convergence_threshold, ...
            'state_file',    state_file, ...
            'save_interval', save_interval, ...
            'resume',        resume_mode);

    detail_matrix = [repmat(m, M_sites, 1), ...
                     repmat(x, M_sites, 1), ...
                     site_indices(:), eta(:)];
    eta_details{idx} = detail_matrix;

    % Trajectories converge at different steps, so align the series to a
    % common time grid before stacking and averaging.
    pop_s_all = align_along_time(population_s_cell);
    pop_s_avg = mean(pop_s_all, 3);
    pop_t_all = align_along_time(totalpopulation_cell);
    pop_t_avg = mean(pop_t_all, 3);

    data(idx, :) = [m, x, mean(eta), t_elapsed];
    results_raw(idx, :) = {population_s_cell, totalpopulation_cell};
    results_avg(idx, :) = {pop_s_avg, pop_t_avg};

    fprintf('Combo m=%d, x=%d done in %.1f s\n', m, x, t_elapsed);
end
toc(tic_all);

all_details = vertcat(eta_details{:});
eta_table = array2table(all_details, ...
    'VariableNames', {'m', 'x', 'site', 'eta'});
disp(eta_table);

figure;
plot(x_r, data(:,3), 'o-', 'LineWidth', 1.5, 'Color', 'b');
xlabel('x'); ylabel('\eta (site-averaged)');
title(sprintf('m = %d,  T_{max} = %.0e', m_f, T_max));
grid on;

figure;
NH = m_f * (m_f+1) .* x_r / 2;
plot(NH, data(:,3), 's-', 'LineWidth', 1.5, 'Color', 'r');
xlabel('NH (number of physical sites)'); ylabel('\eta (site-averaged)');
title(sprintf('m = %d,  T_{max} = %.0e', m_f, T_max));
grid on;

delete(gcp('nocreate'));

function A = align_along_time(cell_arr)
% Stack time series with different row counts into a 3-D array along the
% 3rd dimension; shorter series are padded to the longest length with their
% final (converged) values.
    M    = numel(cell_arr);
    ncol = size(cell_arr{1}, 2);
    nrow = max(cellfun(@(v) size(v, 1), cell_arr));
    A = zeros(nrow, ncol, M, 'like', cell_arr{1});
    for k = 1:M
        v = cell_arr{k};
        t = size(v, 1);
        if t == 0
            continue;
        end
        A(1:t, :, k) = v;
        if t < nrow
            A(t+1:end, :, k) = repmat(v(end, :), nrow - t, 1);
        end
    end
end
