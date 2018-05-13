function u = EVOL_LBF(u0,Img,K,KI,KONE,nu,timestep,mu,lambda1,lambda2,epsilon,numIter)
%  u = EVOL_LBF(u0,Img,Ksigma,KI,KONE,nu,timestep,mu,lambda1,lambda2,epsilon,numIter)

u=u0;
for k1=1:numIter
    u=NeumannBoundCond(u);
    C=curvature_central(u);    % div()  
    HeavU=Heaviside(u,epsilon);
    DiracU=Dirac(u,epsilon);
    
    [f1,f2]=LBF_LocalBinaryFit(K,Img,KI,KONE,HeavU);    
    LBF=LBF_dataForce(Img,K,KONE,f1,f2,lambda1,lambda2);
  
    areaTerm=-DiracU.*LBF;
    penalizeTerm=mu*(4*del2(u)-C);
    lengthTerm=nu.*DiracU.*C;
    u=u+timestep*(lengthTerm+penalizeTerm+areaTerm);
end

% Make a function satisfy Neumann boundary condition
function g = NeumannBoundCond(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

function k = curvature_central(u)
% compute curvature for u with central difference scheme
[ux,uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);
Nx = ux./normDu;
Ny = uy./normDu;
[nxx,junk] = gradient(Nx);
[junk,nyy] = gradient(Ny);
k = nxx+nyy;

function [f1,f2] = LBF_LocalBinaryFit(K,Img,KI,KONE,H)
I=Img.*H;
c1=conv2(H,K,'same');
c2=conv2(I,K,'same');
f1=c2./(c1);
f2=(KI-c2)./(KONE-c1);

function h = Heaviside(x,epsilon)     % function (11)
h=0.5*(1+(2/pi)*atan(x./epsilon));

function f = Dirac(x, epsilon)    % function (12)
f=(epsilon/pi)./(epsilon^2.+x.^2);

function f=LBF_dataForce(Img,K,KONE,f1,f2,lamda1,lamda2)
s1=lamda1.*f1.^2-lamda2.*f2.^2;
s2=lamda1.*f1-lamda2.*f2;
f=(lamda1-lamda2)*KONE.*Img.*Img+conv2(s1,K,'same')-2.*Img.*conv2(s2,K,'same');
