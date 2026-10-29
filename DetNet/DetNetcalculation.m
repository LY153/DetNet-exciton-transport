function [eta, t1, population_s_cell, totalpopulation_cell, c] = ...
    DetNetcalculation(m, x, site_indices, T_max_max, convergence_threshold, ...
                                                 varargin)
% DetNetcalculation  Kinetic Monte-Carlo propagation of the exciton dynamics.
% Runs one trajectory per entry of site_indices and returns the trapping
% efficiency eta, the elapsed time and the recorded population time series.

    p = inputParser;
    addParameter(p, 'state_file', '', @ischar);
    addParameter(p, 'resume', false, @islogical);
    addParameter(p, 'record_path', tempdir, @ischar);
    addParameter(p, 'save_interval', 200, @isscalar);
    addParameter(p, 'buffer_cap', 100, @isscalar);
    parse(p, varargin{:});
    state_file  = p.Results.state_file;
    do_resume   = p.Results.resume;
    rec_path    = p.Results.record_path;
    save_int    = p.Results.save_interval;
    buf_cap     = p.Results.buffer_cap;

    function fname = unique_matfile_path(base_dir, sim)
        rng('shuffle');
        suffix = sprintf('%d_%d_%d', feature('getpid'), randi(1e9), sim);
        fname = fullfile(base_dir, ['sim_', suffix, '.mat']);
    end

    if do_resume && exist(state_file, 'file')
        load(state_file, 'NH', 'NH_full', 'vec_exciton', 'N', 'dt', ...
             'Rpc', 'Rnc', 'U', 'M', 'site_indices', ...
             'converged', 'step_idx', 'Na', 'State', 'c', ...
             'record_files', 'total_steps_done', 'buf_cap');
        record_mats = cell(1, M);
        for sim = 1:M
            record_mats{sim} = matfile(record_files{sim}, 'Writable', true);
        end
        buf_pop = cell(1, M);
        buf_tot = cell(1, M);
        buf_count = 0;
        latest_tot = zeros(M, NH_full, 'single');
        for sim = 1:M
            tmp = record_mats{sim}.totalpopulation;
            if ~isempty(tmp)
                latest_tot(sim, :) = tmp(end, :);
            else
                latest_tot(sim, :) = single(Na(1:NH_full, sim)') / single(N);
            end
        end
        skip_init = true;
    else
        skip_init = false;
        [H, NH] = DetNetGeneration(m, x);
        unitf = 2*pi*3e-5;
        dt = 1;
        N = 1e4;
        N = NH * round(N / NH);

        RT = 1000;
        R_cr = 1e-4;
        R_fl = 1e-6;

        [H02, pR, c_PD, pRedg] = LMR_Para2(H, RT, dt);
        [vec, Hp] = eig(H02);
        vec_d = vec';
        c = vec_d(:, 1)';
        [~, idx] = max(abs(c));
        if c(idx) < 0, c = -c; end
        vec_exciton = vec;

        Hs = (Hp - min(diag(Hp)) * eye(NH)) * unitf;
        t0 = dt:dt:RT;
        Gamma = PureDephasingRate2(c_PD, pRedg, NH, t0);
        disR = ppval(pR, t0);
        clear pR c_PD pRedg H02 Hp vec_d

        NH_full = NH + 2;
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

        LTR = disR(:, :, RT) + diag(Gamma(:, RT));
        clear disR

        R = LTR;
        R_ext = zeros(NH_full);
        R_ext(1:NH, 1:NH) = R;
        R = R_ext;
        clear R_ext LTR
        R(NH_full - 1, 1:NH) = repmat(R_fl, 1, NH);
        R(NH_full, 1:NH) = R_cr * c;

        Rpc = R .* (R > 0);
        Rnc = -R .* (R < 0);

        totR = diag(sum(R, 1));
        HH = Hs - 1i * totR / 2;
        U = expm(-1i * HH * dt);
        clear totR HH

        M = length(site_indices);
        InitStates = zeros(NH_full, M);
        for sim = 1:M
            InitStates(site_indices(sim), sim) = 1;
        end
        State = vec' * InitStates;
        State = State ./ sqrt(sum(abs(State).^2, 1));
        clear InitStates

        Na = zeros(NH_full + 1, M);
        Na(NH_full + 1, :) = N;

        converged = false(M, 1);
        step_idx = 1;
        total_steps_done = 0;

        record_files = cell(1, M);
        record_mats = cell(1, M);
        for sim = 1:M
            fname = unique_matfile_path(rec_path, sim);
            record_files{sim} = fname;
            tmp_mat = matfile(fname, 'Writable', true);
            tmp_mat.population_s = zeros(0, NH, 'single');
            tmp_mat.totalpopulation = zeros(0, NH_full, 'single');
            record_mats{sim} = tmp_mat;
        end

        buf_pop = cell(1, M);
        buf_tot = cell(1, M);
        for sim = 1:M
            buf_pop{sim} = zeros(buf_cap, NH, 'single');
            buf_tot{sim} = zeros(buf_cap, NH_full, 'single');
        end
        buf_count = 0;

        latest_tot = zeros(M, NH_full, 'single');
    end

    if ~skip_init
        for sim = 1:M
            rho_exc0 = Na(NH_full + 1, sim) * (State(1:NH, sim) * State(1:NH, sim)');
            rho_site0 = vec_exciton * rho_exc0 * vec_exciton';
            row_pop = single(real(diag(rho_site0)')) / single(N);
            row_tot = single(Na(1:NH_full, sim)') / single(N);
            buf_count = buf_count + 1;
            buf_pop{sim}(buf_count, :) = row_pop;
            buf_tot{sim}(buf_count, :) = row_tot;
            latest_tot(sim, :) = row_tot;
        end
        total_steps_done = 1;
        if buf_count >= buf_cap
            flush_buffers();
        end
    end

    function flush_buffers()
        if buf_count == 0, return; end
        for s = 1:M
            record_mats{s}.population_s(end+1:end+buf_count, :) = buf_pop{s}(1:buf_count, :);
            record_mats{s}.totalpopulation(end+1:end+buf_count, :) = buf_tot{s}(1:buf_count, :);
        end
        buf_count = 0;
    end

    function save_state_atomically()
        if isempty(state_file), return; end
        flush_buffers();
        tmp_file = [state_file, '.tmp_', sprintf('%d', feature('getpid'))];
        try
            save(tmp_file, 'NH', 'NH_full', 'vec_exciton', 'N', 'dt', ...
                 'Rpc', 'Rnc', 'U', 'M', 'site_indices', ...
                 'converged', 'step_idx', 'Na', 'State', 'c', ...
                 'record_files', 'total_steps_done', 'buf_cap');
            movefile(tmp_file, state_file, 'f');
        catch
            warning('Failed to save state file.');
        end
    end

    tic;
    while step_idx <= T_max_max && ~all(converged)
        active = ~converged;
        Na_old = Na;
        tmpNa = Na_old;

        if any(active)
            State(:, active) = U * State(:, active);
            norms = sqrt(sum(abs(State(:, active)).^2, 1));
            State(:, active) = State(:, active) ./ norms;
        end

        if any(active)
            Pop_active = abs(State(:, active)).^2;
            JumpProb_active = ((R .* (R > 0)) * Pop_active) * dt;
            P_total_active = sum(JumpProb_active, 1);
            n_q_active = Na_old(NH_full + 1, active);
            n_jump_active = binornd(n_q_active, min(P_total_active, 1));

            active_idx = find(active);
            n_active = length(active_idx);
            counts_sub = zeros(NH_full, n_active);

            valid = (P_total_active > 0) & (n_jump_active > 0);
            if any(valid)
                prob_norm = JumpProb_active(:, valid)';
                prob_norm = prob_norm ./ sum(prob_norm, 2);
                counts_valid = mnrnd(n_jump_active(valid)', prob_norm);
                valid_idx_in_active = find(valid);
                counts_sub(:, valid_idx_in_active) = counts_valid';
            end

            tmpNa(1:NH_full, active_idx) = tmpNa(1:NH_full, active_idx) + counts_sub;
            tmpNa(NH_full + 1, active_idx) = tmpNa(NH_full + 1, active_idx) - sum(counts_sub, 1);
            clear counts_sub counts_valid prob_norm
        end

        for sim = 1:M
            if converged(sim), continue; end

            for k = NH_full:-1:1
                nk = Na_old(k, sim);
                if nk <= 0, continue; end
                nk_int = round(nk);
                if nk_int == 0, continue; end

                dPpc = Rpc(:, k) * dt;
                dPnc1 = (Rnc(k, 1:NH_full))' .* (Na_old(1:NH_full, sim) / Na_old(k, sim)) * dt;
                dPnc0 = (Rnc(k, 1:NH_full))' .* (abs(State(1:NH_full, sim)).^2) ...
                        * (Na_old(NH_full + 1, sim) / Na_old(k, sim)) * dt;

                p_jump = [dPpc; dPnc1; dPnc0];
                p_sum = sum(p_jump);
                if p_sum > 1
                    p_all = p_jump / p_sum;
                    no_jump_prob = 0;
                else
                    p_all = p_jump;
                    no_jump_prob = 1 - p_sum;
                end
                p_vec = [p_all; no_jump_prob];

                target_idx = [(1:NH_full)'; (1:NH_full)'; repmat(NH_full + 1, NH_full, 1); k];

                counts = mnrnd(nk_int, p_vec(:)');
                add_counts = accumarray(target_idx, counts(:), [NH_full + 1, 1]);
                tmpNa(1:NH_full + 1, sim) = tmpNa(1:NH_full + 1, sim) + add_counts;
                tmpNa(k, sim) = tmpNa(k, sim) - nk_int;
            end
        end

        Na(:, active) = tmpNa(:, active);

        buf_count = buf_count + 1;
        for sim = 1:M
            row_tot = single(Na(1:NH_full, sim)') / single(N);
            rho_exc = Na(NH_full + 1, sim) * (State(1:NH, sim) * State(1:NH, sim)') ...
                    + diag(Na(1:NH, sim));
            rho_site = vec_exciton * rho_exc * vec_exciton';
            row_pop = single(real(diag(rho_site)')) / single(N);

            buf_pop{sim}(buf_count, :) = row_pop;
            buf_tot{sim}(buf_count, :) = row_tot;
            latest_tot(sim, :) = row_tot;
        end
        total_steps_done = total_steps_done + 1;

        if buf_count >= buf_cap
            flush_buffers();
        end

        for sim = 1:M
            if converged(sim), continue; end
            loss_sum = latest_tot(sim, NH_full - 1) + latest_tot(sim, NH_full);
            if loss_sum >= convergence_threshold
                converged(sim) = true;
            end
        end

        if ~isempty(state_file) && mod(step_idx, save_int) == 0
            flush_buffers();
            save_state_atomically();
        end

        step_idx = step_idx + 1;
    end

    flush_buffers();
    if ~isempty(state_file)
        save_state_atomically();
    end

    t1 = toc;

    population_s_cell = cell(1, M);
    totalpopulation_cell = cell(1, M);
    for sim = 1:M
        population_s_cell{sim} = record_mats{sim}.population_s;
        totalpopulation_cell{sim} = record_mats{sim}.totalpopulation;
    end

    eta = zeros(1, M);
    for sim = 1:M
        eta(sim) = double(totalpopulation_cell{sim}(end, NH_full));
    end
end
