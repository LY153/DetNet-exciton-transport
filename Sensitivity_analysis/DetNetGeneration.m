function [H,NH]=DetNetGeneration(m,x,Jc,energy_gradient_total, ...
                                 dipole_mode,dipole_sigma,energy_disorder)
% DetNetGeneration  Build the network Hamiltonian.
% m is the number of layers (excluding the terminal site); x is the number
% of sites in the second layer. dipole_mode selects fixed parallel dipoles
% ('parallel') or gaussian angular disorder with std dipole_sigma in degrees
% ('gaussian_angle'). Returns H in cm^-1 and the number of sites NH.

% Default values (keep the old call signature working)
if nargin < 3 || isempty(Jc), Jc = 1; end
if nargin < 4 || isempty(energy_gradient_total), energy_gradient_total = 700; end
if nargin < 5 || isempty(dipole_mode), dipole_mode = 'parallel'; end
if nargin < 6 || isempty(dipole_sigma), dipole_sigma = 0; end
if nargin < 7 || isempty(energy_disorder), energy_disorder = 0; end

NH = m*(m+1)*x/2;
Xlen = 124.23; Ylen = 31.97; Zlen = 94.29; % extent of the antenna site box, in units of 10^-1 nm
SiteEnergy = zeros(NH,1); % site energies, in units of cm^-1

% With energy_disorder > 0, each layer draws a uniform offset in
% +/- energy_disorder (cm^-1); all sites of a layer share the same offset.
if energy_disorder > 0
    dE_layer = energy_disorder * (2*rand(m,1) - 1);
else
    dE_layer = zeros(m,1);
end

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

     SiteEnergy(n1:n2,1) = (i-1)*energy_gradient_total/(m-1) + 14100 ...
                           + dE_layer(i);
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

switch dipole_mode
    case 'parallel'
        dipR = repmat([0,1,0],NH,1);
    case 'gaussian_angle'
        % Polar angle gaussian-distributed around [0,1,0] with std
        % dipole_sigma (degrees); azimuth uniform in [0, 2*pi].
        sigma_rad = dipole_sigma * pi / 180;
        theta = abs( sigma_rad * randn(NH,1) );
        theta = min(theta, pi);
        phi = 2*pi*rand(NH,1);
        dipR = [sin(theta).*cos(phi), cos(theta), sin(theta).*sin(phi)];
        dipR = dipR ./ sqrt(sum(dipR.^2,2));
    otherwise
        error('Unknown dipole_mode: %s', dipole_mode);
end

dipR=dipR./repmat(sqrt(sum(dipR.^2,2)),1,3);

D= 4.48; % magnitude of the electric dipole, in units of Debye

H= zeros(NH,NH);
for j= 1: NH-1
    r = repmat(posR(j,:),NH-j,1)-posR(j+1:NH,:);
    normR = sqrt(sum(r.^2,2));
    d1=repmat(dipR(j,:),NH-j,1);
    d2=dipR(j+1:NH,:);
    J = 5042*(D^2)*(dot(d1,d2,2)./(normR.^3)-3*dot(d1,r,2).*dot(d2,r,2)./(normR.^5)); % point-dipole coupling, in units of cm^-1
    H(j, j+1: NH)=J;
end
H = H + H';
H = H + diag(SiteEnergy);
coup = H(2:end,1);
H(1) = 14100-3*max(coup);

H(abs(H)<Jc)=0;
