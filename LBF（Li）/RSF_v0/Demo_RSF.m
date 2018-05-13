% This Matlab file demomstrates a level set method in Chunming Li et al's paper
%    "Minimization of Region-Scalable Fitting Energy for Image Segmentation",
%    IEEE Trans. Image Processing, vol. 17 (10), pp.1940-1949, 2008.
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://www.engr.uconn.edu/~cmli/
%
% Note 1: This method with a small scale parameter sigma, such as sigma = 3, is sensitive to 
%         initialization of the level set function. Appropriate initial level set functions are given in 
%         this code for different test images.
% Note 2: There are several ways to improve the original LBF model to make it robust to initialization.
%         One of the improved LBF algorithms is implemented by the code in the folder LBF_v0.1


clc;clear all;close all;
c0 = 2;
imgID = 1; % 1,2,3,4,5  % choose one of the five test images

Img = imread('1.bmp');%[num2str(imgID),'.bmp']
Img = double(Img(:,:,1));

switch imgID
    case 1
        iterNum =300;
        lambda1 = 1.0;
        lambda2 = 2.0;
        nu = 0.004*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(20:70,30:90) = -c0;
    case 2
        iterNum =200;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.002*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(26:32,28:34) = -c0;
    case 3
        iterNum =500;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.003*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(5:78,32:55) = -c0;
    case 4
        iterNum =150;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.001*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(53:77,46:70) = -c0;
    case 5
        iterNum =220;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.001*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(47:60,86:99) = -c0;
end

u = initialLSF;
figure;imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
title('Initial contour');
[c,h] = contour(u,[0 0],'r');
pause(0.1);

timestep = .1;% time step
mu = 1;% coefficient of the level set (distance) regularization term P(\phi)

epsilon = 1.0;% the papramater in the definition of smoothed Dirac function
sigma=3.0;    % scale parameter in Gaussian kernel
              % Note: A larger scale parameter sigma, such as sigma=10, would make the LBF algorithm more robust 
              %       to initialization, but the segmentation result may not be as accurate as using
              %       a small sigma when there is severe intensity inhomogeneity in the image. If the intensity
              %       inhomogeneity is not severe, a relatively larger sigma can be used to increase the robustness of the LBF
              %       algorithm.
K=fspecial('gaussian',round(2*sigma)*2+1,sigma);     % the Gaussian kernel
I = Img;
KI=conv2(Img,K,'same');     % compute the convolution of the image with the Gaussian kernel outside the iteration
                            % See Section IV-A in the above IEEE TIP paper for implementation.
                                                 
KONE=conv2(ones(size(Img)),K,'same');  % compute the convolution of Gassian kernel and constant 1 outside the iteration
                                       % See Section IV-A in the above IEEE TIP paper for implementation.

% start level set evolution
for n=1:iterNum
    u=RSF(u,I,K,KI,KONE, nu,timestep,mu,lambda1,lambda2,epsilon,1);
    if mod(n,20)==0
        pause(0.1);
        imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
        [c,h] = contour(u,[0 0],'r');
        iterNum=[num2str(n), ' iterations'];
        title(iterNum);
        hold off;
    end
end
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r');
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);

figure;
mesh(u);
title('Final level set function');

