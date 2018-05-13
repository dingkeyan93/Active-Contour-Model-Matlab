function phi = evolution_LGIF(Img,K,KI,KONE,phi,M,lambda1,lambda2,mu,lambda,delta_t,epsilon,numIter)
%   evolution_withoutedge(I, phi0, mu, nu, lambda_1, lambda_2, delta_t, delta_h, epsilon, numIter);
%   input: 
%       I: input image
%       phi0: level set function to be updated
%       mu: weight for length term
%       nu: weight for area term, default value 0
%       lambda_1:  weight for c1 fitting term
%       lambda_2:  weight for c2 fitting term
%       delta_t: time step
%       epsilon: parameter for computing smooth Heaviside and dirac function
%       numIter: number of iterations
%   output: 
%       phi: updated level set function
%  

for k=1:numIter  
    phi=NeumannBoundCond(phi);
    delta_h=Delta(phi,epsilon);
    H=Heaviside(phi,epsilon);
    Curv=curvature(phi);
    [C1,C2]=binaryfit(phi,Img,H,epsilon);
    [f1,f2]=LBF_LocalBinaryFit(K,Img,KI,KONE,H);    
     LBF=LBF_dataForce(Img,K,KONE,f1,f2,lambda1,lambda2);
    % updating the phi function
     F1=(1-M)*LBF;
     F2=M*(-lambda1*(Img-C1).^2+lambda2*(Img-C2).^2);
     penalizingTerm=mu.*(4*del2(phi)-Curv);
     weightedLengthTerm=lambda.*delta_h.*Curv ;
    phi=phi+delta_t*(delta_h.*(F2-F1)+ weightedLengthTerm+penalizingTerm);
end



function g=NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);

function H = Heaviside(phi,epsilon) 
%   Heaviside(phi,epsilon)  compute the smooth Heaviside function
H = 0.5*(1+ (2/pi)*atan(phi./epsilon));

function [f1,f2] = LBF_LocalBinaryFit(K,Img,KI,KONE,H)
I=Img.*H;
c1=conv2(H,K,'same');
c2=conv2(I,K,'same');
f1=c2./(c1);
f2=(KI-c2)./(KONE-c1);

function f=LBF_dataForce(Img,K,KONE,f1,f2,lamda1,lamda2)
s1=lamda1.*f1.^2-lamda2.*f2.^2;
s2=lamda1.*f1-lamda2.*f2;
f=(lamda1-lamda2)*KONE.*Img.*Img+conv2(s1,K,'same')-2.*Img.*conv2(s2,K,'same');

function [C1,C2]= binaryfit(phi,U,H,epsilon) 
a= H.*U;
numer_1=sum(a(:)); 
denom_1=sum(H(:));
C1 = numer_1/denom_1;

b=(1-H).*U;
numer_2=sum(b(:));
c=1-H;
denom_2=sum(c(:));
C2 = numer_2/denom_2;
