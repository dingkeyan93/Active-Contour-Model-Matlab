% The Matlab code of the paper 
% "K. Ding, A Simple Method to improve Initialization Robustness for Active Contours driven by Local Region Fitting Energy, 2018, arXiv:1802.10437."
% Keyan Ding, 2018. 03. 15

clc; clear all; close all;
Img = imread('1.bmp');
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
nu = 0.001*255*255; % the length term. 
lambda1 = 1; 
lambda2 = 1; 
epsilon = 1.0;
timestep = 0.05;
iterNum = 200; % the number of iterations. 
sigma=3; % control the local size

Ksigma=fspecial('gaussian',round(2*sigma)*2+1,sigma);
KONE = imfilter(ones(size(Img)),Ksigma,'replicate');
KI = imfilter(Img,Ksigma,'replicate');
KI2 = imfilter(Img.^2,Ksigma,'replicate');

% --- model weight ---
RSF = 1; 
LRCV = 0;
LIF = 0;
LGD = 0;
CV = 0;

% See "K. Ding, A Simple Method to improve Initialization Robustness for Active Contours driven by Local Region Fitting Energy, 2018, arXiv:1802.10437."
isExchange = 1; % '1' for bright object and dark backgroud; 
                % '-1' for dark object and bright backgroud;
                % '0' represent original model.

% -----start level set evolution-----
h2=figure(2);
tic
for n=1:iterNum
    [u,f1,f2]=ACM(u,Img,Ksigma,KI,KI2,KONE,nu,timestep,mu,epsilon,lambda1,lambda2,CV,RSF,LRCV,LIF,LGD,isExchange);
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

