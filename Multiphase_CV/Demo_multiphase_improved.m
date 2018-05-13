%   This code implements the improved Vese-Chan multiphase level set
%   model in [1] with modification by adding the distance regularizing (DR) term
%   introduced in Li et al's paper [2].
%
%   Reference:
%   [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.
%   [2] Level Set Evolution Without Reinitialization: A New Variational Formulation", CVPR 2005.
%
%   Although reinitialization in Vese and Chan's original level set methods
%   is not required for all images, the level set function can still be
%   degraded after a certain number of iterations. This not only causes the curve irregular, but also slow down 
%   the curve evolution.By adding the distance regularizing term (the internal energy)in Li's level set method [2], 
%   the need for reinitialization is completely eliminated, and the level set function is smooth during the 
%   evolution, which ensure more accurate computation.
%
%   Note that the distance regularizing term has an effect of maintaining the level set function as an approximate 
%   signed distance function near the zero level set. 
%
%   Note: There may be more sophiscated numerical schemes with better performance than the one in
%   this implementation. We only use a simple difference scheme in this version to demonstrate the desirable 
%   distance regularizing effect.
%   Author: Chunming Li

clear;
delta_t = .1;
lambda_1=1;
lambda_2=1;
h = 1;
epsilon=1;
nu = .001*255*255;
fun_n=2;  % two level set functions for 4 phase level set formulation

Img=imread('fourblock_gray.bmp');  
U=Img(:,:,1);

% get the size
[nrow, ncol] =size(U);

ic=nrow/2;
jc=ncol/2;
r=30;
phi = initial_sdf2circle(nrow,ncol,ic,jc,r,fun_n); 
I=double(U);

figure;
imagesc(uint8(I));colormap(gray)
hold on;
plotLevelSet(phi(:,:,1),0,'r');
plotLevelSet(phi(:,:,2),0,'b');

numIter = 10;

for k=1:70
    mu=0.5;  % coefficient for the distance regularizing term (internal energy) in Li's CVPR05 paper
    phi=EVOLUTION_4PHASE_DR(I, phi, nu, lambda_1, lambda_2, mu, delta_t, epsilon, numIter);   % update level set function
    if mod(k,2)==0
        pause(.1);
        imagesc(uint8(I));colormap(gray)
        hold on;
          phi_1=phi(:,:,1);
          phi_2=phi(:,:,2);
          plotLevelSet(phi_1,0,'r');     
          plotLevelSet(phi_2,0,'b'); 
          hold off;
   end           
end
figure;mesh(phi_1);  
title('\phi_1, improved Vese-Chan model by Li');
figure;mesh(phi_2);
title('\phi_2, improved Vese-Chan model by Li');