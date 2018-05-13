function phi = EVOLUTION_CV(I, phi0, mu, nu, lambda_1, lambda_2, delta_t, epsilon, numIter);
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
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li


I=BoundMirrorExpand(I);
phi=BoundMirrorExpand(phi0);

for k=1:numIter
    phi=BoundMirrorEnsure(phi);
    delta_h=Delta(phi,epsilon);
    Curv = curvature(phi);
    [C1,C2]=binaryfit(phi,I,epsilon);
    % updating the phi function
    phi=phi+delta_t*delta_h.*(mu*Curv-nu-lambda_1*(I-C1).^2+lambda_2*(I-C2).^2);    
end
phi=BoundMirrorShrink(phi);