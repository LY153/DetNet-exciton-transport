function [eta, t1, population_s_cell, totalpopulation_cell, c] = ...
    DetNetcalculation(m, x, site_indices, T_max_max, convergence_threshold, ...
                                                 varargin)
% DetNetcalculation  Kinetic Monte-Carlo propagation of the exciton dynamics.
% The independent trajectories (one per entry of site_indices) are run in
% parallel; returns the trapping efficiency eta and the recorded time series.

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

    if do_resume && exist(state_file, 'file')
        load(state_file, 'NH', 'NH_full', 'vec_exciton', 'N', 'dt', ...
             'U', 'Rpc_dt', 'RncT_dt', 'S2', 'target_idx_base', ...
             'M', 'site_indices', 'c', 'record_files', ...
             'Na', 'State', 'converged', 'steps_done');
    else
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

        Gamma_RT = PureDephasingRate2(c_PD, pRedg, NH, RT);
        disR_RT  = ppval(pR, RT);
        disR_RT  = reshape(disR_RT, NH, NH, []);
        LTR = disR_RT(:, :, 1) + diag(Gamma_RT);
        clear pR c_PD pRedg H02 Hp vec_d Gamma_RT disR_RT

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
        clear totR HH R

        Rpc_dt  = Rpc * dt;
        RncT_dt = Rnc.' * dt;
        S2      = abs(vec_exciton).^2;
        target_idx_base = [(1:NH_full)'; (1:NH_full)'; repmat(NH_full+1, NH_full, 1)];
        clear Rpc Rnc

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

        converged  = false(1, M);
        steps_done = zeros(1, M);

        record_files = cell(1, M);
        for sim = 1:M
            fname = unique_matfile_path(rec_path, sim);
            record_files{sim} = fname;
            tmp_mat = matfile(fname, 'Writable', true);
            tmp_mat.population_s    = zeros(0, NH, 'single');
            tmp_mat.totalpopulation = zeros(0, NH_full, 'single');
        end
    end

    tic;
    Na_out        = zeros(NH_full + 1, M);
    State_out     = zeros(NH_full, M);
    converged_out = false(1, M);
    steps_out     = zeros(1, M);

    parfor sim = 1:M
        [conv_s, Na_s, St_s, st_s] = run_one_sim( ...
            record_files{sim}, converged(sim), ...
            Na(:, sim), State(:, sim), steps_done(sim), ...
            T_max_max, convergence_threshold, NH, NH_full, N, dt, ...
            U, Rpc_dt, RncT_dt, vec_exciton, S2, target_idx_base, ...
            save_int, buf_cap);
        converged_out(sim) = conv_s;
        Na_out(:, sim)     = Na_s;
        State_out(:, sim)  = St_s;
        steps_out(sim)     = st_s;
    end
    t1 = toc;

    Na = Na_out; State = State_out; converged = converged_out; steps_done = steps_out;

    population_s_cell    = cell(1, M);
    totalpopulation_cell = cell(1, M);
    for sim = 1:M
        mf = matfile(record_files{sim});
        population_s_cell{sim}    = mf.population_s;
        totalpopulation_cell{sim} = mf.totalpopulation;
    end

    eta = zeros(1, M);
    for sim = 1:M
        eta(sim) = double(totalpopulation_cell{sim}(end, NH_full));
    end

    if ~isempty(state_file)
        try
            tmp_file = [state_file, '.tmp_', sprintf('%d', feature('getpid'))];
            save(tmp_file, 'NH', 'NH_full', 'vec_exciton', 'N', 'dt', ...
                 'U', 'Rpc_dt', 'RncT_dt', 'S2', 'target_idx_base', ...
                 'M', 'site_indices', 'c', 'record_files', ...
                 'Na', 'State', 'converged', 'steps_done');
            movefile(tmp_file, state_file, 'f');
        catch
            warning('Failed to save state file.');
        end
    end
end

function fname = unique_matfile_path(base_dir, sim)
    rng('shuffle');
    suffix = sprintf('%d_%d_%d', feature('getpid'), randi(1e9), sim);
    fname = fullfile(base_dir, ['sim_', suffix, '.mat']);
end

function [converged_sim, Na_sim, State_sim, steps_done] = run_one_sim( ...
        rec_file, converged_in, Na_in, State_in, steps_done_in, ...
        T_max_max, convergence_threshold, NH, NH_full, N, dt, ...
        U, Rpc_dt, RncT_dt, vec_exciton, S2, target_idx_base, ...
        save_int, buf_cap)
% run_one_sim - run a single Monte-Carlo trajectory (one initial site).
% All inputs are passed by value, so this local function is safe to call
% from inside a parfor loop.

    converged_sim = converged_in;
    Na_sim        = Na_in;
    State_sim     = State_in;
    steps_done    = steps_done_in;

    if converged_sim
        return;
    end

    mf = matfile(rec_file, 'Writable', true);

    buf_pop   = zeros(buf_cap, NH, 'single');
    buf_tot   = zeros(buf_cap, NH_full, 'single');
    buf_count = 0;

    if steps_done == 0
        w0 = vec_exciton * State_sim(1:NH);
        row_pop = single(real(Na_sim(NH_full+1)*(w0.*conj(w0)) + S2*Na_sim(1:NH))') / single(N);
        row_tot = single(Na_sim(1:NH_full)') / single(N);
        buf_count = 1;
        buf_pop(1, :) = row_pop;
        buf_tot(1, :) = row_tot;
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

        row_tot = single(Na_sim(1:NH_full)') / single(N);
        w = vec_exciton * State_sim(1:NH);
        row_pop = single(real(Na_sim(NH_full+1)*(w.*conj(w)) + S2*Na_sim(1:NH))') / single(N);

        buf_count = buf_count + 1;
        buf_pop(buf_count, :) = row_pop;
        buf_tot(buf_count, :) = row_tot;
        steps_done = steps_done + 1;

        if buf_count >= buf_cap
            mf.population_s(end+1:end+buf_count, :)    = buf_pop(1:buf_count, :);
            mf.totalpopulation(end+1:end+buf_count, :) = buf_tot(1:buf_count, :);
            buf_count = 0;
        end

        loss_sum = row_tot(NH_full - 1) + row_tot(NH_full);
        if loss_sum >= convergence_threshold
            converged_sim = true;
        end

        if mod(steps_done, save_int) == 0 && buf_count > 0
            mf.population_s(end+1:end+buf_count, :)    = buf_pop(1:buf_count, :);
            mf.totalpopulation(end+1:end+buf_count, :) = buf_tot(1:buf_count, :);
            buf_count = 0;
        end
    end

    if buf_count > 0
        mf.population_s(end+1:end+buf_count, :)    = buf_pop(1:buf_count, :);
        mf.totalpopulation(end+1:end+buf_count, :) = buf_tot(1:buf_count, :);
    end
end
