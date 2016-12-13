function Zdot = myodefunc(t,Z,L,A,n)
% fast system (*)\prime
% Z(1)=phi/sigma, Z(2)=n*v/sigma, Z(3)=u/phi
% phi = \theta^((1+a)/(1+n))
Zdot = zeros(3,1);
p = Z(1);
q = Z(2);
r = Z(3);

D = 1+2*A-n;
c0 = 2 / D;
d1 = -2*(A-n) / D;
b  = 1/D + (1+n)/D*L;

Zdot(1) = n^1 * ( p * ( (1+A)/(1+n)/L*( r^(1+n)-c0 ) - (d1+L*p*r+q) ) );
Zdot(2) = n^1 * ( q * ( 1-q-L*p*r ) + b*p*r );
Zdot(3) = n^0 * ( r * ( (A-n)/(1+n)/L*( r^(1+n)-c0 ) + (d1+L*p*r+q) ) );

end