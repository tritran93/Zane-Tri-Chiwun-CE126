function [phi_prime] = double_pendulum(t,phi)

phi_prime = zeros(4,1);

%Data setup

m1=2;
m2=1;
L1=1;
L2=2;
g = 9.81;
F1 = 0;
F2 = 0;
f1 = F1*cos(omega1*t);
f2 = F2*cos(omega2*t);

%Assign names for coefficient
%Reminder phi(1) = phi1; phi(2) = phi1.dot;
%phi(3) = phi2; phi(4) = phi2.dot;
a = (m1+m2)*L1;                             
b = m2*L2*cos(phi(1)-phi(3));               
c = - m2*L2*phi(4)*phi(4)*sin(phi(1)-phi(3)) - (m1+m2)*g*sin(phi(1))+f1;     
d = m2*L1*cos(phi(1)-phi(3));
e = m2*L2;
f = m2*L1*sin(phi(1)-phi(3)) - m2*g*sin(phi(3)) + f2;

%Setting up the function based on the matrix. See homework note for work.
phi_prime(1) = phi(2);
phi_prime(3) = phi(4);
phi_prime(2) = (e*c - b*f)/(a*e - b*d);
phi_prime(4) = (f*a - d*c)/(e*a - b*d);

end
