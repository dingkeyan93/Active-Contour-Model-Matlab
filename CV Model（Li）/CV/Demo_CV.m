%   Matlad code implementing Chan-Vese model in the paper 'Active Contours Without Edges'
%   Author: Chunming Li, all right reserved
%   email: li_chunming@hotmail.com
%   URL:   http://www.engr.uconn.edu/~cmli

clear all;
close all;
Img=imread('three.bmp');   
% Img=imread('vessel3.bmp'); % Note: this an example of images with intensity inhomogeneity. 
                             % CV model does not work for this image. 
Img=double(Img(:,:,1));

% get the size
[nrow,ncol] =size(Img);

ic=nrow/2;
jc=ncol/2;
r=20;
initialLSF = sdf2circle(nrow,ncol,ic,jc,r);
u=initialLSF;


numIter = 150;
timestep = 0.1;
lambda_1=1;
lambda_2=1;
h = 1;
epsilon=1;
nu = 0.001*255*255;  % tune this parameter for different images

figure;
imagesc(Img,[0 255]);colormap(gray)
hold on;
contour(u,[0 0],'r');


% start level set evolution
for k=1:numIter
    u=EVOL_CV(Img, u, nu, lambda_1, lambda_2, timestep, epsilon, 1);   % update level set function
    if mod(k,10)==0
        pause(.1);
        imagesc(Img,[0 255]);colormap(gray)
        hold on;
        contour(u,[0 0],'r');
        hold off;
    end    
end;
