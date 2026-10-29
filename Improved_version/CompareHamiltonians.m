function result = CompareHamiltonians(m, x, varargin)
% CompareHamiltonians  Compare the eigenvalues and eigenstates of the
% Hamiltonian before and after the coupling cutoff.
%
% H0 is the Hamiltonian before H(abs(H)<Jc)=0 is applied in
% DetNetGeneration, H is the one after it. Both are rebuilt here
% independently, with exactly the same rules as DetNetGeneration, and the
% original file is not modified.
%
% Input:
%   m, x - same meaning as in DetNetGeneration.
% Optional name-value parameters:
%   'overlap_tol' - threshold on |<v0_i, v_i>|; indices below it are
%                   flagged for re-matching (default 0.9).
%   'n_search'    - number of nearby eigenstates searched when
%                   re-matching (default 10).
%   'Jc'          - cutoff coupling constant in cm^-1 (default 1).
%   'verbose'     - print the comparison table (default true).
% Output (result struct):
%   E0, E     - eigenvalues of H0 and H (sorted in descending order).
%   dE        - E - E0 (eigenvalue shift caused by the cutoff).
%   overlap   - |<v0_i, v_i>| for identically sorted indices.
%   dev       - 1 - overlap.
%   bad_idx   - indices flagged for re-matching.
%   corrected - re-matching table (if any).
%   NH, Jc    - number of sites, cutoff coupling constant.
%
% Example:
%   result = CompareHamiltonians(5, 6);

    p = inputParser;
    addParameter(p, 'overlap_tol', 0.9, @(v) isscalar(v) && v > 0 && v < 1);
    addParameter(p, 'n_search',    10,  @(v) isscalar(v) && v > 0);
    addParameter(p, 'Jc',           1,  @(v) isscalar(v) && v > 0);
    addParameter(p, 'verbose',   true,  @islogical);
    parse(p, varargin{:});
    overlap_tol = p.Results.overlap_tol;
    n_search    = p.Results.n_search;
    Jc          = p.Results.Jc;
    verbose     = p.Results.verbose;

    [H0, H, Jc] = build_hamiltonians(m, x, Jc);
    NH = size(H0, 1);

    % Consistency check: DetNetGeneration uses the fixed Jc = 1, so H must
    % match its output only when Jc == 1.
    if Jc == 1
        H_ref = DetNetGeneration(m, x);
        assert(max(abs(H(:) - H_ref(:))) < 1e-9, ...
            'Constructed H differs from DetNetGeneration; check build_hamiltonians.');
    end

    [V0, E0] = eig(H0);  E0 = diag(E0);
    [V, E] = eig(H);  E = diag(E);
    [E0, idx0] = sort(E0, 'descend');  V0 = V0(:, idx0);
    [E, idx] = sort(E, 'descend');  V = V(:, idx );

    dE = E - E0;

    overlap = abs(diag(V0' * V)).^2;
    dev     = 1 - overlap;

    bad  = find(overlap < overlap_tol);
    corrected = zeros(0, 6);
    for ii = 1:numel(bad)
        i = bad(ii);
        [~, order] = sort(abs(E - E0(i)));
        cand = order(1:min(n_search, numel(order)));
        ov   = abs(V0(:, i)' * V(:, cand)).^2;
        [ov_max, k] = max(ov);
        j_best = cand(k);
        corrected = [corrected; i, E0(i), j_best, E(j_best), ov_max, 1 - ov_max]; %#ok<AGROW>
    end

    % Best match within an energy window for every eigenstate: this removes
    % energy-order swaps caused by the cutoff and the arbitrary basis choice
    % inside degenerate subspaces, so dev_corrected reflects only the genuine
    % mixing of each eigenstate caused by the cutoff.
    overlap_corrected = zeros(NH, 1);
    for i = 1:NH
        [~, order] = sort(abs(E - E0(i)));
        cand = order(1:min(n_search, numel(order)));
        overlap_corrected(i) = max(abs(V0(:, i)' * V(:, cand)).^2);
    end
    dev_corrected = 1 - overlap_corrected;

    if verbose
        fprintf('\n========== m = %d, x = %d, NH = %d, Jc = %g ==========\n', m, x, NH, Jc);
        fprintf('Eigenvalues sorted in descending order (H0: before cutoff, H: after cutoff):\n');
        fprintf('%6s %16s %16s %16s %12s %12s\n', ...
            'No.', 'E0 (H0)', 'E (H)', 'dE = E-E0', 'overlap', '1-overlap');
        for i = 1:NH
            fprintf('%6d %16.6f %16.6f %16.6f %12.6f %12.6f\n', ...
                i, E0(i), E(i), dE(i), overlap(i), dev(i));
        end

        if isempty(bad)
            fprintf('\nAll matched-state overlaps are >= %.3f; no re-matching needed.\n', overlap_tol);
        else
            fprintf('\nEntries with |<v0_i, v_i>| < %.3f were re-matched among nearby eigenstates:\n', overlap_tol);
            fprintf('%6s %16s %10s %16s %12s %12s\n', ...
                'i', 'E0(i)', 'j_best', 'E(j_best)', 'ov_corr', '1-ov_corr');
            for r = 1:size(corrected, 1)
                fprintf('%6d %16.6f %10d %16.6f %12.6f %12.6f\n', ...
                    corrected(r,1), corrected(r,2), corrected(r,3), ...
                    corrected(r,4), corrected(r,5), corrected(r,6));
            end
            fprintf('If ov_corr is much larger than the original overlap, the cutoff\n');
            fprintf('swapped the energy ordering of the eigenstates; if it is still far\n');
            fprintf('below 1, the eigenstate itself is strongly mixed by the cutoff.\n');
        end
    end

    result.E0      = E0;
    result.E       = E;
    result.dE      = dE;
    result.overlap = overlap;
    result.dev     = dev;
    result.overlap_corrected = overlap_corrected;
    result.dev_corrected     = dev_corrected;
    result.NH      = NH;
    result.Jc      = Jc;
    result.bad_idx = bad;
    if isempty(corrected)
        result.corrected = [];
    else
        result.corrected = array2table(corrected, ...
            'VariableNames', {'i', 'E0_i', 'j_best', 'E_jbest', 'overlap_corr', 'dev_corr'});
    end
end

function [H0, H, Jc] = build_hamiltonians(m, x, Jc)
% Same construction as DetNetGeneration, but also returns H0 (before the
% cutoff). Jc is the cutoff coupling constant in cm^-1.
    if nargin < 3 || isempty(Jc)
        Jc = 1;
    end
    NH = m*(m+1)*x/2;
    Xlen = 124.23; Ylen = 31.97; Zlen = 94.29;
    posR = zeros(NH,3);
    SiteEnergy = zeros(NH,1);

    for i = 1:m
        n = i*x; n1 = 1+i*(i-1)*x/2; n2 = i*(i+1)*x/2;
        t_fine = linspace(0, 2*pi, 10000);
        ds_dt = sqrt((i*(Xlen/(2*m))*cos(t_fine)).^2 + (i*(Zlen/(2*m))*sin(t_fine)).^2);
        s_fine = cumtrapz(t_fine, ds_dt);
        L = s_fine(end);
        s_target = (0:n-1) * L / n;
        theta = interp1(s_fine, t_fine, s_target);
        posR(n1:n2,2) = i*Ylen/m;
        posR(n1:n2,1) = i*(Xlen/(2*m))*sin(theta);
        posR(n1:n2,3) = i*(Zlen/(2*m))*cos(theta);
        SiteEnergy(n1:n2,1) = (i-1)*700/(m-1) + 14100;
    end
    NH = NH + 1;
    posR1 = zeros(NH,3);
    posR1(1,:) = [0,0,0];
    posR1(2:end,:) = posR;
    posR = posR1;
    SiteEnergy1 = zeros(NH,1);
    SiteEnergy1(1,1) = 14000;
    SiteEnergy1(2:end,1) = SiteEnergy;
    SiteEnergy = SiteEnergy1;
    dipR = repmat([0,1,0],NH,1);
    dipR = dipR./repmat(sqrt(sum(dipR.^2,2)),1,3);
    D = 4.48;

    H0 = zeros(NH,NH);
    for j = 1:NH-1
        r = repmat(posR(j,:),NH-j,1)-posR(j+1:NH,:);
        normR = sqrt(sum(r.^2,2));
        d1 = repmat(dipR(j,:),NH-j,1);
        d2 = dipR(j+1:NH,:);
        J = 5042*(D^2)*(dot(d1,d2,2)./(normR.^3)-3*dot(d1,r,2).*dot(d2,r,2)./(normR.^5));
        H0(j, j+1:NH) = J;
    end
    H0 = H0 + H0';
    H0 = H0 + diag(SiteEnergy);
    coup = H0(2:end,1);
    H0(1) = 14100-3*max(coup);
    H = H0;
    H(abs(H) < Jc) = 0;
end
