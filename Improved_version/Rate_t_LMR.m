function [pR,cRpd,Rpd]=Rate_t_LMR(E, coe_mat, Reg, Redg, Reddg, Img, Imdg, Imddg, dtg, dtt)
% Rate_t_LMR  Time-dependent modified-Redfield energy-transfer rates.
% INPUT:
%   E       - exciton energy vector, in fs^-1
%   coe_mat - eigenvector matrix of the system Hamiltonian (columns are
%             eigenstates). The overlap coefficients
%                 J(p,q,r,s) = sum_n conj(C(n,p))*C(n,q)*conj(C(n,r))*C(n,s)
%             are computed on the fly, avoiding the huge 4D Jcoe array.
%   Reg / Img     - lineshape function, real / imaginary part
%   Redg / Imdg   - first derivative of the lineshape, real / imaginary part
%   Reddg / Imddg - time correlation function, real / imaginary part
%   dtg, dtt      - time steps of the arrays above
% OUTPUT:
%   pR   - rate matrix pR(b,a,t) in fs^-1, splined over time; rates run from
%          the column index to the row index, e.g. pR(1,2) is the rate from
%          exciton 2 to exciton 1
%   cRpd - coefficient matrix of the time-dependent pure-dephasing rate,
%          Rpd(t) = cRpd * Redg(t)
%   Rpd  - pure-dephasing rate matrix evaluated at the final bath time

    n = size(coe_mat,1);
    tmax=(size(Reg,2)-1)*dtg;

    tt = 0:dtt:tmax;

    nt=size(tt,2);

    dG =Redg+1i*Imdg;
    L  =-Imdg(:,end); % in fs^-1

    dR = zeros(n,n,nt);
    cRpd = zeros(n,n);
    Rpd= zeros(n,n);

    C  = coe_mat;
    cc = conj(C);

    for a = 1:n
        ca  = C(:,a);
        cca = cc(:,a);
        for b = [1:a-1 a+1:n]
            cb  = C(:,b);
            ccb = cc(:,b);

            J_bbaa = sum(ccb.*cb.*cca.*ca);
            J_bbbb = sum(ccb.*cb.*ccb.*cb);
            J_aaaa = sum(cca.*ca.*cca.*ca);
            J_abbb = sum(cca.*cb.*ccb.*cb);
            J_abaa = sum(cca.*cb.*cca.*ca);
            J_babb = sum(ccb.*ca.*ccb.*cb);
            J_baaa = sum(ccb.*ca.*cca.*ca);
            J_abba = sum(cca.*cb.*ccb.*ca);

            c  = conj(2*J_bbaa - J_bbbb - J_aaaa);
            va = conj(J_abbb - J_abaa);
            vb = conj(J_babb - J_baaa);

            V = conj(J_abba)*(Reddg+1i*Imddg)...
                    -(va*vb)*dG.^2 ...
                    +1i*2*((va*conj(J_baaa)+vb*conj(J_abaa)).*L)*dG...
                    +4*(conj(J_baaa)*conj(J_abaa))*(L.^2);

            dR(b,a,:)=2*real(exp(1i*(E(a)-E(b))*tt+c*Reg(1:dtt/dtg:end)+...
                1i*c*(L*tt+Img(1:dtt/dtg:end))).*V);
            cRpd(b,a) = -c;
            Rpd(b,a)= -c*Redg(:,end);

        end
    end

    pR=spline(tt, dR);
    clear dR

    pR.order = 5;
    coe = pR.coefs;
    ndt = pR.pieces;

    pR.coefs = [coe(:,1)/4 coe(:,2)/3 coe(:,3)/2 coe(:,4) zeros(size(coe,1), 1)];
    coe = pR.coefs;
    coefs_sum = [coe(:,1)*dtt^4 coe(:,2)*dtt^3 coe(:,3)*dtt^2 coe(:,4)*dtt];
    pR.coefs = reshape(pR.coefs, n*n, ndt, 5);
    coefs_sum = reshape(coefs_sum, n*n, ndt, 4);
    pR.coefs(:, 2:ndt,5)= cumsum(sum(coefs_sum(:, 1:ndt-1, :), 3),2);
    pR.coefs = reshape(pR.coefs, n*n*ndt, 5);

end
