function phi = EVOLUTION_4PHASE(I, phi0, nu, lambda_1, lambda_2, delta_t, epsilon, numIter);
%   This function updates the level set function in the Vese-Chan multiphase level set
%   model in [1].
%
%   Reference:
%   [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.
%  
%   created on 01/20/2006 by Chunming Li and Lu Li
%   author: Chunming Li and Lu Li
%   email: li_chunming@hotmail.com
%   Homepage: 
%            http://www.vuiis.vanderbilt.edu/~licm/
%            http://www.engr.uconn.edu/~cmli
%   Copyright (c) 2004-2008 by Chunming Li and Lu Li

phi(:,:,1)=phi0(:,:,1);
phi(:,:,2)=phi0(:,:,2);

for m=1:numIter
    for k=1:2
        tmp_phi=phi(:,:,k);
        tmp_phi=NeumannBoundCond(tmp_phi);
        delta_h=Delta(tmp_phi,epsilon);  
        differenceScheme=0; % use 0, 1, 2 for different schemes to compute curvature
        Curv = CURVATURE_CV(tmp_phi,differenceScheme); 
        [C,mult]=quadrifit(phi,I,epsilon,2);
        b=zeros(size(Curv));

        if k==1
            b=(C(1)-C(2))*(C(1)+C(2)-2*I).*mult(:,:,1,1);
            b=b+(C(3)-C(4))*(C(3)+C(4)-2*I).*mult(:,:,1,2);
        else
            b=(C(1)-C(3))*(C(1)+C(3)-2*I).*mult(:,:,2,1);
            b=b+(C(2)-C(4))*(C(2)+C(4)-2*I).*mult(:,:,2,2);            
        end
        tmp_phi=tmp_phi+delta_t*(delta_h.*(nu*Curv+b));    
        phi(:,:,k)=tmp_phi;        
    end
end

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);          