function [eta, t1, c] = ...
    DetNetcalculation(m, x, site_indices, T_max_max, convergence_threshold, ...
                                                 varargin)
% DetNetcalculation  Kinetic Monte-Carlo propagation with disorder and
% ensemble splitting: N_r Hamiltonian realizations are built and the N-system
% ensemble is split over them. Returns the ensemble-averaged efficiency eta
% (one value per initial site) and the elapsed time.

    p = inputParser;
    addParameter(p, 'state_file', '', @ischar);
    addParameter(p, 'resume', false, @islogical);
    addParameter(p, 'R_cr', 1e-4, @isscalar);
    addParameter(p, 'R_fl', 1e-6, @isscalar);
    addParameter(p, 'Jc', 1, @isscalar);
    addParameter(p, 'energy_gradient_total', 700, @isscalar);
    addParameter(p, 'dipole_mode', 'parallel', @ischar);
    addParameter(p, 'dipole_sigma', 0, @isscalar);
    addParameter(p, 'initial_state_mode', 'site_average', @ischar);
    addParameter(p, 'N_r', 1, @isscalar);
    addParameter(p, 'energy_disorder', 0, @isscalar);
    addParameter(p, 'rng_seed', [], ...
                 @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v >= 0));
    parse(p, varargin{:});
    state_file  = p.Results.state_file;
    do_resume   = p.Results.resume;
    R_cr = p.Results.R_cr;
    R_fl = p.Results.R_fl;
    Jc = p.Results.Jc;
    energy_gradient_total = p.Results.energy_gradient_total;
    dipole_mode = p.Results.dipole_mode;
    dipole_sigma = p.Results.dipole_sigma;
    initial_state_mode = p.Results.initial_state_mode;
    N_r = p.Results.N_r;
    energy_disorder = p.Results.energy_disorder;
    rng_seed = p.Results.rng_seed;
    if ~(isnumeric(N_r) && isscalar(N_r) && N_r >= 1 && mod(N_r, 1) == 0)
        error('N_r must be a positive integer.');
    end
    if ~(isnumeric(energy_disorder) && isscalar(energy_disorder) && ...
         energy_disorder >= 0)
        error('energy_disorder must be a non-negative scalar.');
    end

    can_resume = do_resume && N_r == 1 && exist(state_file, 'file');
    if can_resume
        try
            load(state_file, 'NH', 'NH_full', 'dt', 'N', 'N_part', 'N_r', ...
                 'M', 'site_indices', 'target_idx_base', ...
                 'U_cell', 'Rpc_dt_cell', 'RncT_dt_cell', ...
                 'c', 'Na', 'State', 'converged', 'steps_done', ...
                 'R_cr', 'R_fl', 'Jc', 'energy_gradient_total', ...
                 'dipole_mode', 'dipole_sigma', ...
                 'initial_state_mode', 'energy_disorder', 'rng_seed');
        catch
            warning('DetNetcalculation:resume', ...
                    ['Existing state file is incompatible with the current ' ...
                     'code; restarting from scratch.']);
            can_resume = false;
        end
        if can_resume && ~isequal(energy_disorder, p.Results.energy_disorder)
            warning('DetNetcalculation:resume', ...
                    ['Resumed state was built with energy_disorder = %g; ' ...
                     'the requested value %g is ignored.'], ...
                    energy_disorder, p.Results.energy_disorder);
        end
    end

    if can_resume
        Na_cell = cell(1, N_r);
        State_cell = cell(1, N_r);
        for r = 1:N_r
            cols = (r - 1) * M + (1:M);
            Na_cell{r}    = Na(:, cols);
            State_cell{r} = State(:, cols);
        end
        rng_base_seed = rng_seed;   % informational only when resuming
        if isempty(rng_base_seed)   % MC streams still need a valid base seed
            rng('shuffle');
            rng_base_seed = randi(2^31 - 1);
        end
    else
        if do_resume && N_r > 1
            warning('DetNetcalculation:resume', ...
                    'Resume is disabled when N_r > 1; restarting from scratch.');
        end

        unitf = 2*pi*3e-5;
        dt = 1;
        N = 2e4;

        RT = 1000;

        % One base seed; part r uses seed (rng_base_seed + r - 1), so every
        % part draws independent random numbers for the random parts of the
        % Hamiltonian (per-layer energy disorder and gaussian_angle dipoles),
        % and the run is reproducible when rng_seed is given.
        if isempty(rng_seed)
            rng('shuffle');
            rng_base_seed = randi(2^31 - 1);
        else
            rng_base_seed = double(rng_seed);
        end

        rng(rng_base_seed);
        [H0, NH] = DetNetGeneration(m, x, Jc, energy_gradient_total, ...
                                    dipole_mode, dipole_sigma, ...
                                    energy_disorder);
        N = NH * round(N / NH);
        if N_r > N
            error('N_r (%d) must not exceed N (%d).', N_r, N);
        end
        N_part = floor(N / N_r) * ones(1, N_r);
        N_part(N_r) = N - sum(N_part(1:N_r-1));

        NH_full = NH + 2;
        target_idx_base = [(1:NH_full)'; (1:NH_full)'; repmat(NH_full+1, NH_full, 1)];

        if strcmp(initial_state_mode, 'site_average')
            M = length(site_indices);
            InitStates = zeros(NH_full, M);
            for sim = 1:M
                InitStates(site_indices(sim), sim) = 1;
            end
        elseif strcmp(initial_state_mode, 'equal_superposition')
            M = 1;
            InitStates = zeros(NH_full, M);
            InitStates(site_indices, 1) = 1 / sqrt(length(site_indices));
        else
            error('Unknown initial_state_mode: %s', initial_state_mode);
        end

        U_cell       = cell(1, N_r);
        Rpc_dt_cell  = cell(1, N_r);
        RncT_dt_cell = cell(1, N_r);
        State_cell   = cell(1, N_r);
        Na_cell      = cell(1, N_r);

        for r = 1:N_r
            if r == 1
                H = H0;
            else
                rng(rng_base_seed + r - 1);
                [H, ~] = DetNetGeneration(m, x, Jc, energy_gradient_total, ...
                                          dipole_mode, dipole_sigma, ...
                                          energy_disorder);
            end

            [H02, pR, c_PD, pRedg] = LMR_Para2(H, RT, dt);
            [vec, Hp] = eig(H02);
            vec_d = vec';
            c_part = vec_d(:, 1)';
            [~, idx] = max(abs(c_part));
            if c_part(idx) < 0, c_part = -c_part; end

            Hs = (Hp - min(diag(Hp)) * eye(NH)) * unitf;

            Gamma_RT = PureDephasingRate2(c_PD, pRedg, NH, RT);
            disR_RT  = ppval(pR, RT);
            disR_RT  = reshape(disR_RT, NH, NH, []);
            LTR = disR_RT(:, :, 1) + diag(Gamma_RT);
            clear pR c_PD pRedg H02 Hp vec_d Gamma_RT disR_RT

            Hs_ext = zeros(NH_full);
            Hs_ext(1:NH, 1:NH) = Hs;
            Hs = Hs_ext;
            clear Hs_ext
            vec_ext = zeros(NH_full);
            vec_ext(1:NH, 1:NH) = vec;
            vec = vec_ext;
            vec(NH_full - 1, NH_full - 1) = 1;
            vec(NH_full, NH_full) = 1;
            clear vec_ext

            R = LTR;
            R_ext = zeros(NH_full);
            R_ext(1:NH, 1:NH) = R;
            R = R_ext;
            clear R_ext LTR
            R(NH_full - 1, 1:NH) = repmat(R_fl, 1, NH);
            R(NH_full, 1:NH) = R_cr * c_part;

            Rpc = R .* (R > 0);
            Rnc = -R .* (R < 0);

            totR = diag(sum(R, 1));
            HH = Hs - 1i * totR / 2;
            U = expm(-1i * HH * dt);
            clear totR HH R

            Rpc_dt  = Rpc * dt;
            RncT_dt = Rnc.' * dt;
            clear Rpc Rnc

            State = vec' * InitStates;
            State = State ./ sqrt(sum(abs(State).^2, 1));

            Na = zeros(NH_full + 1, M);
            Na(NH_full + 1, :) = N_part(r);

            U_cell{r}       = U;
            Rpc_dt_cell{r}  = Rpc_dt;
            RncT_dt_cell{r} = RncT_dt;
            State_cell{r}   = State;
            Na_cell{r}      = Na;
        end
        c = c_part;  % returned for compatibility; only meaningful when N_r == 1

        converged  = false(1, N_r * M);
        steps_done = zeros(1, N_r * M);
    end

    tic;
    n_tasks       = N_r * M;
    Na_out        = zeros(NH_full + 1, n_tasks);
    State_out     = zeros(NH_full, n_tasks);
    converged_out = false(1, n_tasks);
    steps_out     = zeros(1, n_tasks);

    parfor task = 1:n_tasks
        % Per-iteration seeding makes the Monte-Carlo noise reproducible
        % (parfor workers otherwise use unrelated global streams).
        rng(rng_base_seed + 1e6 + task);
        r   = ceil(task / M);
        sim = task - (r - 1) * M;
        [conv_s, Na_s, St_s, st_s] = run_one_sim( ...
            converged(task), ...
            Na_cell{r}(:, sim), State_cell{r}(:, sim), steps_done(task), ...
            T_max_max, convergence_threshold, NH_full, N_part(r), ...
            U_cell{r}, Rpc_dt_cell{r}, RncT_dt_cell{r}, target_idx_base);
        converged_out(task) = conv_s;
        Na_out(:, task)     = Na_s;
        State_out(:, task)  = St_s;
        steps_out(task)     = st_s;
    end
    t1 = toc;

    Na = Na_out; State = State_out; converged = converged_out; steps_done = steps_out;

    % Efficiency: Na_out(NH_full, task) is the number of systems captured by
    % the reaction center for part r and initial site sim. Summing over the
    % N_r parts and dividing by the total ensemble size N gives the
    % system-count-weighted average over the whole ensemble.
    eta = zeros(1, M);
    for r = 1:N_r
        cols = (r - 1) * M + (1:M);
        eta = eta + Na_out(NH_full, cols);
    end
    eta = eta / N;

    if ~isempty(state_file) && N_r == 1
        try
            tmp_file = [state_file, '.tmp_', sprintf('%d', feature('getpid'))];
            save(tmp_file, 'NH', 'NH_full', 'dt', 'N', 'N_part', 'N_r', ...
                 'M', 'site_indices', 'target_idx_base', ...
                 'U_cell', 'Rpc_dt_cell', 'RncT_dt_cell', ...
                 'c', 'Na', 'State', 'converged', 'steps_done', ...
                 'R_cr', 'R_fl', 'Jc', 'energy_gradient_total', ...
                 'dipole_mode', 'dipole_sigma', ...
                 'initial_state_mode', 'energy_disorder', 'rng_seed', ...
                 'rng_base_seed');
            movefile(tmp_file, state_file, 'f');
        catch
            warning('Failed to save state file.');
        end
    end
end

function [converged_sim, Na_sim, State_sim, steps_done] = run_one_sim( ...
        converged_in, Na_in, State_in, steps_done_in, ...
        T_max_max, convergence_threshold, NH_full, N, ...
        U, Rpc_dt, RncT_dt, target_idx_base)
% run_one_sim - run a single Monte-Carlo trajectory (one part x initial site).
% Only the final reaction-center population (efficiency) is needed, so no
% per-step time series is recorded here.
% All inputs are passed by value, so this local function is safe to call
% from inside a parfor loop.

    converged_sim = converged_in;
    Na_sim        = Na_in;
    State_sim     = State_in;
    steps_done    = steps_done_in;

    if converged_sim
        return;
    end

    while steps_done < T_max_max && ~converged_sim
        State_sim = U * State_sim;
        State_sim = State_sim ./ sqrt(sum(abs(State_sim).^2, 1));

        Na_old = Na_sim;
        tmpNa  = Na_old;

        Pop      = abs(State_sim).^2;
        JumpProb = Rpc_dt * Pop;
        P_total  = sum(JumpProb);
        n_q      = Na_old(NH_full + 1);
        n_jump   = binornd(n_q, min(P_total, 1));
        if n_jump > 0 && P_total > 0
            prob_norm = JumpProb' / P_total;
            counts    = mnrnd(n_jump, prob_norm);
            tmpNa(1:NH_full)   = tmpNa(1:NH_full) + counts';
            tmpNa(NH_full + 1) = tmpNa(NH_full + 1) - n_jump;
        end

        for k = NH_full:-1:1
            nk = Na_old(k);
            if nk <= 0, continue; end
            nk_int = round(nk);
            if nk_int == 0, continue; end

            dPpc  = Rpc_dt(:, k);
            dPnc1 = RncT_dt(:, k) .* (Na_old(1:NH_full) / nk);
            dPnc0 = RncT_dt(:, k) .* Pop * (Na_old(NH_full + 1) / nk);

            p_jump = [dPpc; dPnc1; dPnc0];
            p_sum  = sum(p_jump);
            if p_sum > 1
                p_all = p_jump / p_sum;
                no_jump_prob = 0;
            else
                p_all = p_jump;
                no_jump_prob = 1 - p_sum;
            end
            p_vec = [p_all; no_jump_prob];

            counts     = mnrnd(nk_int, p_vec(:)');
            add_counts = accumarray([target_idx_base; k], counts(:), [NH_full + 1, 1]);
            tmpNa      = tmpNa + add_counts;
            tmpNa(k)   = tmpNa(k) - nk_int;
        end

        Na_sim = tmpNa;

        steps_done = steps_done + 1;

        loss_sum = (Na_sim(NH_full - 1) + Na_sim(NH_full)) / N;
        if loss_sum >= convergence_threshold
            converged_sim = true;
        end
        if mod(steps_done, 1e4) == 0
            fprintf('  step=%d, loss_sum=%.4f, eta_now=%.4f\n', ...
                    steps_done, loss_sum, Na_sim(NH_full) / N);
        end
    end
end
