function C=ell_general(kc,p,c,s)

Np=length(kc);
C=zeros(1,Np);

for ii=1:Np
    
fun = @(x) (c*cos(x).^2+s*sin(x).^2)./((cos(x).^2+p*sin(x).^2).*sqrt(cos(x).^2+kc(ii)^2*sin(x).^2)) ;
    C(ii) = integral(fun,0,pi/2);
end

% fun = @(x) (c*cos(x).^2+s*sin(x).^2)./((cos(x).^2+p*sin(x).^2).*sqrt(cos(x).^2+kc.^2*sin(x).^2)) ;
%     C= integral(fun,0,pi/2);

end


