%%% For Staggered reference and AP layers
%%% Behtash Behin-Aein
function [Hcpl] = Hcpl_Disc_Offset(delta,Ms,R_free,t_free,crt)
global R t
% Thickness of the Ref or bottom magnet
% t

%%% tolerance
tol = 1e-6;
fac=1-1e-9;

%%% distance between magnets
% Delta

czt=t_free/2+delta+t/2;

%dblquad
qy = dblquad(@By_of_single_magnet,czt-t_free/2,czt+t_free/2,crt,crt+fac*R_free,tol);
Dyy=qy/(R_free*t_free);

convert=Ms*10^6*(4*pi*1e-7)*10 *1e3;

Hcpl = Dyy*convert;


