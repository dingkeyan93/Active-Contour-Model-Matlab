% This Matlab file demomstrates a level set algorithm based on Chunming Li et al's paper:
% "Implicit Active Contours Driven By Local Binary Fitting Energy" in Proceedings of CVPR'07


clc;clear all;close all;
c0 =2;
imgID=3;

Img=imread('vessel.bmp');
%Img=imread('vessel2.bmp');  % uncommont this line to use ther other vessel image
I=Img(:,:,1);
Img=double(Img);

switch imgID
     case 1
       phi= ones(size(Img(:,:,1))).*c0;
       a=43;b=51;c=20;d=28;
       phi(a:b,c:d) = -c0;
       figure;
       imshow(I);colormap;
       hold on;
       plotLevelSet(phi, 0, 'g');
       hold off;
    case 2
       [m,n]=size(Img(:,:,1));
       a=m/2; b=n/2;r=5;
       phi= ones(m,n).*c0;
       phi(a-r:a+r,b-r:b+r) = -c0;
       imshow(I);colormap;
       hold on;
       plotLevelSet(phi, 0, 'r');
       hold off;
    case 3
       figure;imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;axis equal;
       text(6,6,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'g');
       BW=roipoly;
       phi=c0*2*(0.5-BW);
       hold on;
       [c,h] = contour(phi,[0 0],'r');
       hold off;
end
pause(0.01);

%²ÎÊýÑ¡Ôñ
iterNum = 400;
lambda1 = 1.0;
lambda2 = 1.0;
nu = 0.002*255*255;
timestep = 0.1;
mu = 1;
epsilon = 1.0;


% scale parameter in Gaussian kernel
sigma=3.0;    
K=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
KI=conv2(Img,K,'same');  
KONE=conv2(ones(size(Img)),K,'same');


% start level set evolution
time = cputime;
for n=1:iterNum
   numIter=1;
    %level set evolution.  
    phi=EVOL_LBF(phi,Img,K,KI,KONE,nu,timestep,mu,lambda1,lambda2,epsilon,numIter);
     if mod(n,10)==0
        pause(0.001);
        imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
        contour(phi,[0 0],'r');
        iterNum=[num2str(n), ' iterations'];
        title(iterNum);
        hold off;
    end
end
totaltime = cputime - time
imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
contour(phi,[0 0],'r');
iterNum=[num2str(n), ' iterations'];
title(iterNum);