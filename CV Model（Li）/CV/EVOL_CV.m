function phi = EVOL_CV(I, phi0, nu, lambda_1, lambda_2, timestep, epsilon, numIter);
%   This function updates the level set function according to the CV model 
%   input: 
%       I: input image
%       phi0: level set function to be updated
%       mu: weight for length term
%       nu: weight for area term, default value 0
%       lambda_1:  weight for c1 fitting term
%       lambda_2:  weight for c2 fitting term
%       muP: weight for level set regularization term 
%       timestep: time step
%       epsilon: parameter for computing smooth Heaviside and dirac function
%       numIter: number of iterations
%   output: 
%       phi: updated level set function
%  
%   created on 04/26/2004
%   Author: Chunming Li, all right reserved
%   email: li_chunming@hotmail.com
%   URL:   http://www.engr.uconn.edu/~cmli/research/

phi=phi0;
for k=1:numIter
    phi=NeumannBoundCond(phi);
    diracPhi=Delta(phi,epsilon);
    Hphi=Heaviside(phi, epsilon);
    kappa = CURVATURE(phi,'cc');
    [C1,C2]=binaryfit(I,Hphi);
    % updating the phi function
    phi=phi+timestep*(diracPhi.*(nu-lambda_1*(I-C1).^2+lambda_2*(I-C2).^2));    
end


function H = Heaviside(phi,epsilon) 
H = 0.5*(1+ (2/pi)*atan(phi./epsilon));

function Delta_h = Delta(phi, epsilon)
Delta_h=(epsilon/pi)./(epsilon^2+ phi.^2);

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  