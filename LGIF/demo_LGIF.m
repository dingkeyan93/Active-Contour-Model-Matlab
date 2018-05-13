%   Matlab code implementing Chan-Vese model in the paper 'Active Contours Without Edges'
%   This method works well for bimodal images, for example the image 'three.bmp'


clc;clear all;close all;
c0 =2;
imgID=3;

Img=imread('cq18.jpg');    
U=Img(:,:,1);

% the initial level set
switch imgID
     case 1
       phi= ones(size(Img(:,:,1))).*c0;
       a=43;b=51;c=20;d=28;
       phi(a:b,c:d) = -c0;
       figure;
       imshow(Img);colormap;
       hold on;
       plotLevelSet(phi, 0, 'g');
       hold off;
    case 2
       [m,n]=size(Img(:,:,1));
       a=m/2; b=n/2;r=5;
       phi= ones(m,n).*c0;
       phi(a-r:a+r,b-r:b+r) = -c0;
       imshow(Img);colormap;
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
pause(1);

%参数选择
pause(0.1);
lambda=0.001*255*255; 
lambda1=1.0;
lambda2=1.0;
delta_t =0.1;
epsilon=1;
mu =1;
numIter = 1;

% CV和LBF的权重
Img=double(U);
M=0.05; 

% scale parameter in Gaussian kernel
sigma=10.0;   
K=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
KI=conv2(Img,K,'same');  
KONE=conv2(ones(size(Img)),K,'same');


% start level set evolution
time = cputime;
for k=1:100,
   phi = evolution_LGIF(Img,K,KI,KONE,phi,M,lambda1,lambda2,mu,lambda,delta_t,epsilon,numIter);  % update level set function
    if mod(k,10)==0
        pause(0.01);
        imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
        contour(phi,[0 0],'r');
        iterNum=[num2str(k), ' iterations'];
        title(iterNum);
        hold off;
    end    
end;
totaltime = cputime - time
imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
contour(phi,[0 0],'r');
iterNum=[num2str(k), ' iterations'];
title(iterNum);
