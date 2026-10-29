function [ Gamma ] = PureDephasingRate2( c_PD , pRedg , NH , t )
% PureDephasingRate2  Markovian pure-dephasing rates Gamma at time t.
A = zeros(NH,1);
M = 0.5*ones(NH);
for j = 1 : NH
    A(j,1) = sum(c_PD(j,j+1:NH))+sum(c_PD(1:NH-1,j));

    M(j,j) = 0.5*(2*NH-j);
    M(j,j+1:NH-1) = 1;
end

Gamma = (M\A)*ppval(pRedg,t);

end
