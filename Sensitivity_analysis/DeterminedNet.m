% DeterminedNet  DetNet driver with disorder and ensemble options: runs the
% simulation for m = m_f layers and x = x_r, saves the efficiency results to
% the current folder, and plots eta vs x and vs NH.
% clear; close all;

m_f = 5;
x_r = 4:8;

T_max = 1e6;

convergence_threshold = 0.97;
resume_mode   = false;             % true = resume computation, false = restart computation

% Baseline parameters; edit these to run sensitivity scans.
R_cr = 1e-4;
R_fl = 1e-6;
Jc = 1;
energy_gradient_total = 700;       % total energy gradient across layers (cm^-1)
dipole_mode = 'gaussian_angle';    % 'parallel' or 'gaussian_angle'
dipole_sigma = 1;                  % std of the gaussian angular disorder (degrees), only for dipole_mode='gaussian_angle'
energy_disorder = 0;               % per-layer energy offset drawn uniformly from +/- energy_disorder (cm^-1)
rng_seed = [];                     % [] = random each run; an integer makes the whole run reproducible
initial_state_mode = 'site_average'; % 'site_average' or 'equal_superposition'
N_r = 20;                          % number of ensemble parts (resume is disabled when N_r > 1)


numCombos = length(x_r);

data = zeros(numCombos, 4);
eta_details = cell(numCombos, 1);

delete(gcp('nocreate'));
try
    n_cores = feature('numcores');
catch
    n_cores = 12;
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

    if isempty(rng_seed)
        seed_tag = 'rnd';
    else
        seed_tag = num2str(rng_seed);
    end
    state_file = fullfile(tempdir, sprintf('state_m%d_x%d_%s_dipsigma%s.mat', ...
                         m, x, initial_state_mode, num2str(dipole_sigma)));

    [eta, t_elapsed, ~] = ...
    DetNetcalculation(m, x, site_indices, ...
        T_max, convergence_threshold, ...
        'state_file',    state_file, ...
        'resume',        resume_mode, ...
        'R_cr',          R_cr, ...
        'R_fl',          R_fl, ...
        'Jc',            Jc, ...
        'energy_gradient_total', energy_gradient_total, ...
        'dipole_mode',   dipole_mode, ...
        'dipole_sigma',  dipole_sigma, ...
        'energy_disorder', energy_disorder, ...
        'rng_seed',      rng_seed, ...
        'initial_state_mode', initial_state_mode, ...
        'N_r',            N_r);

    n_eta = length(eta);
    if n_eta == M_sites
        detail_matrix = [repmat(m, M_sites, 1), ...
                         repmat(x, M_sites, 1), ...
                         site_indices(:), eta(:)];
    else
        detail_matrix = [repmat(m, n_eta, 1), ...
                         repmat(x, n_eta, 1), ...
                         (1:n_eta)', eta(:)];
    end
    eta_details{idx} = detail_matrix;

    data(idx, :) = [m, x, mean(eta), t_elapsed];

    fprintf('Combo m=%d, x=%d done in %.1f s\n', m, x, t_elapsed);
end
toc(tic_all);

all_details = vertcat(eta_details{:});
eta_table = array2table(all_details, ...
    'VariableNames', {'m', 'x', 'site', 'eta'});
disp(eta_table);

summary_table = array2table(data, ...
    'VariableNames', {'m', 'x', 'mean_eta', 't_elapsed'});
out_base = fullfile(pwd, sprintf('eta_results_m%d_%s_Nr%d_ED%s', m_f, ...
                    initial_state_mode, N_r, num2str(energy_disorder)));
mat_out = [out_base '.mat'];
sum_csv = [out_base '_summary.csv'];
det_csv = [out_base '_details.csv'];
save(mat_out, 'data', 'summary_table', 'eta_table', 'all_details', ...
     'm_f', 'x_r', 'T_max', 'convergence_threshold', ...
     'R_cr', 'R_fl', 'Jc', 'energy_gradient_total', ...
     'dipole_mode', 'dipole_sigma', ...
     'energy_disorder', 'rng_seed', 'initial_state_mode', 'N_r');
writetable(summary_table, sum_csv);
writetable(eta_table, det_csv);
fprintf('Efficiency results saved to:\n  %s\n  %s\n  %s\n', ...
        mat_out, sum_csv, det_csv);

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
