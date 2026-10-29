function [Reg,Reintcf,Retcf,Img,Imintcf,Imtcf]=bath_construction(J, F, KbT, t_max, dtg, dtt)
% bath_construction  Construct the bath correlation functions and lineshape.
% Outputs: Retcf/Imtcf - real/imaginary parts of the time correlation
%          function; Reintcf/Imintcf - their integrals; Reg/Img - real and
%          imaginary parts of the lineshape function.
% t_max must be long enough for the coupling strength, e.g. t_max > 15000
% for coupling <= 5.

%Time Axis Assignment
tg = 0:dtg:t_max;

tt = 0:dtt:t_max;

%REAL PART
%The first point of the frequency grid is avoided because coth() is singular
%there; the w = 0 contribution is added analytically (tcf_w0).
dww=1/10/t_max*dtt;
ww=(dww:dww:1/dtt)*2*pi;

Redtcf=@(w) J(w).*coth(w/(2*KbT));
tcf_w0=2*F(0)*KbT;

tcf_fft = 2*pi*(1/dtt)*real(ifft([tcf_w0 Redtcf(ww)]));
clear ww Redtcf
Retcf=tcf_fft(1:t_max/dtt+1);
clear tcf_fft

%TCF integral Re. part
ptcf=spline(tt, Retcf);
pintcf = ptcf;
coe = pintcf.coefs;
ndt   = pintcf.pieces;
pintcf.order = 5;

pintcf.coefs = [coe(:,1)/4 coe(:,2)/3 coe(:,3)/2 coe(:,4) zeros(ndt, 1)];
coe = pintcf.coefs;
coefs_sum = [coe(:,1)*dtt^4 coe(:,2)*dtt^3 coe(:,3)*dtt^2 coe(:,4)*dtt];
pintcf.coefs(2:ndt,5)= cumsum(sum(coefs_sum(1:ndt-1, :),2),1);
clear coe coefs_sum

Reintcf=ppval(tt, pintcf);

%Lineshape Function Re. part
gstep =ppval(tg, pintcf)*dtg+0.5*ppval(tg, ptcf)*dtg^2;
Reg = cumsum(gstep, 2);
clear gstep ptcf pintcf

%IMAGINARY PART
Imdtcf =@(w) -J(w);

ww =(0:1/2/t_max:1/dtt)*2*pi;

Imddg_ifft = 2*pi/dtt*ifft(Imdtcf(ww));
Imddg_tt = imag(Imddg_ifft);

Imtcf=Imddg_tt(1:t_max/dtt+1);
clear Imddg_ifft Imdtcf

%TCF integral Im. part
Imdintcf =@(w) F(w);

Imdg_ifft = 2*pi/dtt*ifft(Imdintcf(ww));
Imdg_tt = real(Imdg_ifft);

Imintcf=Imdg_tt(1:t_max/dtt+1)-real(Imdg_ifft(1))*ones(1,t_max/dtt+1);
clear Imdg_ifft Imdintcf ww

%Lineshape Function Im. part
gstep =spline(tt, Imintcf, tg)*dtg+0.5*spline(tt, Imtcf, tg)*dtg^2;
Img = cumsum(gstep, 2);
clear gstep

end
