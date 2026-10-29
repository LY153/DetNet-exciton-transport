function [ H , pR , c_PD , pRedg ] = LMR_Para2(H,t_max,dt)
% LMR_Para2  Prepare the bath and the modified-Redfield rates for the
% Hamiltonian H (in cm^-1). Rate outputs are in fs^-1.
ne = size(H,1);

% Unit conversion factor from cm^-1 to fs^-1
unitf=2*pi*3e-5;

% Bath parameters [wc Lo de dis], in cm^-1:
%   wc  = Ohmic cut-off frequency
%   Lo  = Ohmic reorganization energy
%   de  = FWHM of the site-energy distribution
%   dis = disorder strength (0 = no disorder)
bath = [100 35 0 0];

Kb=0.6950*unitf;
T = 300; % temperature (K)

dtg = dt;
dtt = dt;

wc = bath(1)*unitf;
Lo = bath(2)*unitf;

% Ohmic spectral density
J=@(w) Lo/wc*w.*exp(-w/wc);
F=@(w) Lo/wc*exp(-w/wc);

[Reg,Redg,Reddg ,Img ,Imdg ,Imddg]=bath_construction(J, F, Kb*T, ...
    t_max, dtg, dtt);

% Exciton basis and exciton energies
[coe_mat,Ee]=eig(H);
Ee=diag(Ee);
W=Ee*unitf;

% The 4D overlap tensor Jcoe is not built explicitly; Rate_t_LMR computes
% the required coefficients on the fly from coe_mat.
[pR,c_PD,~] = Rate_t_LMR(W, coe_mat, Reg, Redg, Reddg, Img, Imdg, Imddg, dtg, dtt);

pRedg=spline(0:dtt:t_max, Redg);

end
