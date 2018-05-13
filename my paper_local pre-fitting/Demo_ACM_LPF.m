% The Matlab code of the paper 
% "K. Ding L. Xiao and G. Weng, Active Contours driven by Local Pre-Fitting Energy for Fast Image Segmentation, Pattern Recognition Letters, 104, 29-36, 2017."
% Keyan Ding, 2018.02.12

clc; clear all; close all;
Img = imread('2.bmp');
Img = double(Img(:,:,1));

% -----set initial contour-----
c0 = 1;
initialLSF = ones(size(Img(:,:,1))).*c0;
initialLSF(15:35,40:60) = -c0;

% figure;imagesc(Img, [0, 255]); colormap(gray);hold on; axis off;axis equal;
% text(6,6,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'g');
% BW=roipoly;
% hold off;
% initialLSF=c0*2*(0.5-BW);

% BW=im2bw(imread('mask.bmp'),0.5);
% initialLSF=c0*2*(0.5-BW);

u = initialLSF;
h1 = figure(1);imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
contour(initialLSF,[0 0],'g','linewidth',1.5);hold off
pause(0.1);

% -----set parameters-----
mu = 1; % the distance regularization term
nu = 0.01*255*255; % the length term. 
lambda1 = 1; 
lambda2 = lambda1; 
epsilon = 1.0;
timestep = 0.02;
iterNum = 200; % the number of iterations. 
sigma = 2; % control the local size

% --- Local pre-fitting functions ---    
K=fspecial('gaussian',round(2*sigma)*2+1,sigma);  
[f1,f2,Im]=LPF(Img,K);
e1=Img.*Img.*imfilter(ones(size(Img)),K,'replicate')-2.*Img.*imfilter(f1,K,'replicate')+imfilter(f1.^2,K,'replicate');
e2=Img.*Img.*imfilter(ones(size(Img)),K,'replicate')-2.*Img.*imfilter(f2,K,'replicate')+imfilter(f2.^2,K,'replicate');
% e1=(Img-f1).^2;
% e2=(Img-f2).^2;

%------color image segmentation------
% Img = double(Img);
% for i = 1:3
%     [f1,f2,Im]=LPF(Img(:,:,i),K);
%     e1(:,:,i)=Img(:,:,i).*Img(:,:,i).*imfilter(ones(size(Img(:,:,i))),K,'replicate')-2.*Img(:,:,i).*imfilter(f1,K,'replicate')+imfilter(f1.^2,K,'replicate');
%     e2(:,:,i)=Img(:,:,i).*Img(:,:,i).*imfilter(ones(size(Img(:,:,i))),K,'replicate')-2.*Img(:,:,i).*imfilter(f2,K,'replicate')+imfilter(f2.^2,K,'replicate');
% end
% e1=-(e1(:,:,1)+e1(:,:,2)+e1(:,:,3))./3;
% e2=-(e2(:,:,1)+e2(:,:,2)+e2(:,:,3))./3;

% -----start level set evolution-----
h2=figure(2);
tic
for n=1:iterNum
    u = ACM_LPF(u,nu,timestep,mu,epsilon,lambda1,lambda2,e1,e2);
    if mod(n,10)==0       
        imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
        contour(u,[0 0],'r'); title([num2str(n), ' iterations ']);
        hold off;pause(.01);
    end
end
toc

% -----display result-----
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r','linewidth',1.5);
% figure;mesh(u);colorbar;title('Final level set function');hold on, contour(u,[0 0],'r','linewidth',1.5);

