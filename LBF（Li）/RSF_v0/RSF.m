function u = RSF(u0,Img,Ksigma,KI,KONE,nu,timestep,mu,lambda1,lambda2,epsilon,numIter)
% LSE_LBF implements the level set evolution (LSE) for the method in Chunming Li et al's paper:
%       "Minimization of Region-Scalable Fitting Energy for Image Segmentation", 
%        IEEE Trans. Image Processing, vol. 17 (10), pp.1940-1949, 2008.
%
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://www.engr.uconn.edu/~cmli/
%
% For easy understanding of my code, please read the comments in the code that refer
% to the corresponding equations in the above IEEE TIP paper. 


u=u0;
for k1=1:numIter
    u=NeumannBoundCond(u);
    K=curvature_central(u);                             

    DrcU=(epsilon/pi)./(epsilon^2.+u.^2);               % eq.(9)

    [f1, f2] = localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon);


    %%% compute lambda1*e1-lambda2*e2
    s1=lambda1.*f1.^2-lambda2.*f2.^2;                   % compute lambda1*e1-lambda2*e2 in the 1st term in eq. (15) in IEEE TIP 08
    s2=lambda1.*f1-lambda2.*f2;
    dataForce=(lambda1-lambda2)*KONE.*Img.*Img+conv2(s1,Ksigma,'same')-2.*Img.*conv2(s2,Ksigma,'same');
                                                        % eq.(15)
    A=-DrcU.*dataForce;                                 % 1st term in eq. (15)
    P=mu*(4*del2(u)-K);                                 % 3rd term in eq. (15), where 4*del2(u) computes the laplacian (d^2u/dx^2 + d^2u/dy^2)
    L=nu.*DrcU.*K;                                      % 2nd term in eq. (15)
    u=u+timestep*(L+P+A);                               % eq.(15)
end

function [f1, f2]= localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon)
% compute f1 and f2
Hu=0.5*(1+(2/pi)*atan(u./epsilon));                     % eq.(8)

I=Img.*Hu;
c1=conv2(Hu,Ksigma,'same');                             
c2=conv2(I,Ksigma,'same');                              % the numerator of eq.(14) for i = 1
f1=c2./(c1);                                            % compute f1 according to eq.(14) for i = 1
f2=(KI-c2)./(KONE-c1);                                  % compute f2 according to the formula in Section IV-A, 
                                                        % which is an equivalent expression of eq.(14) for i = 2.
                                                         

function g = NeumannBoundCond(f)
% Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

function k = curvature_central(u)                       
% compute curvature
[ux,uy] = gradient(u);                                  
normDu = sqrt(ux.^2+uy.^2+1e-10);                       % the norm of the gradient plus a small possitive number 
                                                        % to avoid division by zero in the following computation.
Nx = ux./normDu;                                       
Ny = uy./normDu;
[nxx,junk] = gradient(Nx);                              
[junk,nyy] = gradient(Ny);                              
k = nxx+nyy;                                            % compute divergence


