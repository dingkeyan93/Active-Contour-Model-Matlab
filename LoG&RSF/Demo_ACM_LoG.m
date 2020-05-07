% The Matlab code of the paper 
% "K. Ding, L. Xiao and G. Weng. Active contours driven by region-scalable fitting and optimized loglacian of Gaussian energy for image segmentation. Signal Processing, 134, 224-233, 2017.
% Keyan Ding, 2018.02.29

clc; clear all; close all;
Img = imread('3.bmp');
Img = double(Img(:,:,1));

% -----set initial contour-----
c0 = 1;
initialLSF = ones(size(Img(:,:,1))).*c0;
initialLSF(15:55,40:80) = -c0;

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
mu = 2; % the distance regularization term
nu = 0.002*255*255; % the length term. (need adjust)
theta = 20; %the LoG term, + / - (need adjust)
lambda1 = 1; 
lambda2 = 1; 
epsilon = 1.0;
timestep = 0.05;
iterNum = 300; % the number of iterations. 
sigma = 3; % control the local size

Ksigma=fspecial('gaussian',round(2*sigma)*2+1,sigma);
KONE = imfilter(ones(size(Img)),Ksigma,'replicate');
KI = imfilter(Img,Ksigma,'replicate');

% --- model weight ---
RSF = 1; 
LRCV = 0;
LIF = 0;

% see "K. Ding, A Simple Method to improve Initialization Robustness for Active Contours driven by Local Region Fitting Energy, 2018, arXiv:1802.10437."
isExchange = 0; % '1' for bright object and dark backgroud; 
                % '-1' for dark object and bright backgroud;
                % '0' represent original model.

% ----- regularised loglacian-----
G=fspecial('gaussian',[9 9],1);
Img_gao=imfilter(Img,G,'replicate');
[Ix,Iy]=gradient(Img_gao);
f=Ix.^2+Iy.^2;
log=4*del2(Img_gao);
% figure,imshow(log,[]),colorbar
g_=zeros(size(Img));
log_=zeros(size(Img));
g=exp(-0.01*f); % alpha = 0.01
for i=1:100
    log_=log_+0.01*(g.*log_ - (1-g).*(log_ - 5*log)); % beta = 5
end
log_ = theta*log_; 
% figure,imshow(log_,[]);colorbar

% -----start level set evolution-----
h2=figure(2);
tic
for n=1:iterNum
    [u,f1,f2]=ACM_LoG(u,Img,Ksigma,KI,KONE,nu,timestep,mu,epsilon,lambda1,lambda2,RSF,LRCV,LIF,isExchange,log_);
    if mod(n,10)==0       
        imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
        contour(u,[0 0],'r'); title([num2str(n), ' iterations ']);
        hold off;pause(.1);
    end
end
toc

% -----display result-----
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r','linewidth',1.5);
% figure;mesh(u);colorbar;title('Final level set function');hold on, contour(u,[0 0],'r','linewidth',1.5);

