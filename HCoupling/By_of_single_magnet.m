function [Bz]=By_of_single_magnet(y,rho)
global R t
b=t/2;   

yp=y+b;
ym=y-b;

beta_p=yp./sqrt(yp.^2+(rho+R).^2);
beta_m=ym./sqrt(ym.^2+(rho+R).^2);

gamma=(R-rho)./(R+rho);

kp=sqrt((yp.^2+(R-rho).^2)/(yp.^2+(R+rho).^2));
km=sqrt((ym.^2+(R-rho).^2)/(ym.^2+(R+rho).^2));


C3=ell_general(kp,gamma^2,1,gamma);
C4=ell_general(km,gamma^2,1,gamma);


Bz=1/pi* R./(R+rho).* (beta_p*C3-beta_m*C4);

end

