function [phi,f1,f2,Hphi] = LIF_2D(I,phi,timestep,epsilon,K)

phi = NeumannBoundCond(phi);

Hphi = Heaviside(phi,epsilon);
DiracPhi = Delta(phi,epsilon);

[f1,f2] = Local_Avr(I,Hphi,K);

phi = phi + timestep*DiracPhi.*((I - f1.*Hphi - f2.*(1 - Hphi)).*(f1 - f2));

function H = Heaviside(phi,epsilon)
H = 0.5*(1+(2/pi)*atan(phi./epsilon));

function Delta_h = Delta(phi,epsilon)
Delta_h = (epsilon/pi)./(epsilon^2+phi.^2);

function [f1,f2] = Local_Avr(I,H,K)

f1 = conv2(I.*H,K,'same');
c1 = conv2(H,K,'same');
f1 = f1./c1;
f2 = conv2(I.*(1-H),K,'same');
c2 = conv2(1-H,K,'same');
f2 = f2./c2;

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);