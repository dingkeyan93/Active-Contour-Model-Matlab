%   This code implements the Vese-Chan multiphase level set model in [1].
%   Note: There may be more sophiscated numerical schemes with better performance than the one in
%   this implementation.
%
%   Reference:
%   [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.

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

for k=1:120
    phi=EVOLUTION_4PHASE(I, phi, nu, lambda_1, lambda_2, delta_t, epsilon, numIter);   % update level set function
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
title('\phi_1, Vese-Chan model');
figure;mesh(phi_2);
title('\phi_2, Vese-Chan model');