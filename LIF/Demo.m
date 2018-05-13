% This Matlab program demomstrates the level set algorithm in paper:
%    "Active contours driven by local image fitting energy" 
%    to appear in Pattern Recognition, 2010
% Author: Kaihua Zhang, Huihui Song and Lei Zhang
% E-mail: zhkhua@mail.ustc.edu.cn, cslzhang@comp.polyu.edu.hk  
% http://www4.comp.polyu.edu.hk/~cslzhang/

%  Notes:
%   1. Some parameters may need to be modified for different types of images. Please contact the author if any problem regarding the choice of parameters.
%   2. Intial contour should be set properly.

% Date: 5/11/2009

clear all;close all;clc;

Img = imread('1.bmp');
Img = double(Img(:,:,1));
sigma =3;% the key parameter which needs to be tuned properly.
sigma_phi = 0.9;% the variance of regularized Gaussian kernel
K = fspecial('gaussian',2*round(2*sigma)+1,sigma);
K_phi = fspecial('gaussian',5,sigma_phi);

[nrow,ncol] = size(Img);

phi = ones(nrow,ncol);
phi(35:nrow-35,55:ncol-35) = -1;
figure;  imagesc(Img,[0 255]);colormap(gray);hold on;
contour(phi,[0 0],'b');

timestep = 1;
epsilon = 1.5;
for n = 1:600   
      [phi,f1,f2,Hphi]= LIF_2D(Img,phi,timestep,epsilon,K);
      phi = conv2(phi,K_phi,'same');
      if mod(n,40)==0
      pause(0.0001);
      imagesc(Img,[0 255]);colormap(gray)
      hold on;contour(phi,[0 0],'b');
      iterNum=[num2str(n), ' iterations'];        
      title(iterNum);        
      hold off;
      end
end

figure;
mesh(phi);