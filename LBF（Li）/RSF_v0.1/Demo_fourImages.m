% This code demomstrates an improved algorithm based on the local binary fitting (LBF) model
% in Chunming Li et al's paper:
%    "Minimization of Region-Scalable Fitting Energy for Image Segmentation", 
%     IEEE Trans. Image Processing, vol. 17 (10), pp.1940-1949, 2008.
%
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://www.engr.uconn.edu/~cmli/
%
%  Notes:
%    1. Some parameters are set to default values for the demos in this package. They may need to be
%       modified for different types of images. 
%    2. The current version does not work for images with multiple junctions, due to its two-phase
%       formulation (i.e. using only one level set function). For example, an image has 3 objects/regions, 
%       and each object/region is directly contiguous to all the other two objects/regions. This code will be 
%       extended to a multiphase formulation in a new version.
%    3. The image intensities may need to be rescaled to the range of [0, 255], if the intensities are much lower
%       or much higher than 255. Alternatively, you can change the parameters lambda1 and lambda2, and nu (the 
%       coefficient of lenght term) accordingly.


close all;
clear all;

imgID=3;   % choose one of the four images by setting imgID to 1, 2, 3, or 4.           
c0 = 2;
if imgID ==1
    iterNum=120;
    lambda1 = 1.0;  lambda2 = 1.0; nu = 0.001*255*255;
    load Vessel2_initialLSF2.mat;   % load image and initial contour
    u=c0*sign(initialLSF);
elseif imgID == 2
    iterNum=200;
    lambda1 = 1.0;  lambda2 = 1.0; nu = 0.001*255*255;
    load Vessel3_initialLSF;  % load image and initial contour
    u=c0*sign(initialLSF);
elseif imgID == 3
    iterNum =200;
    lambda1 = 1.0;  lambda2 = 2.0; nu = 0.003*255*255;
    load  mybrainC100b_initialLSF.mat; % load image and initial contour
    u=initialLSF;
elseif imgID==4
    iterNum=50;
    lambda1 = 1.0;  lambda2 = 1.0; nu = 0.001*255*255;
    load nonhomo_initialLSF.mat; % load image and initial contour
    u=c0*sign(initialLSF);
end

Img=double(Img(:,:,1));
figure;imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
[nrow, ncol]=size(Img);
contour(u,[0 0],'r');
title('Initial LSF')
pause(.05);


pause(0.1);
timestep = .1;
mu = 1;

epsilon = 1.0;
sigma=3.0;    % scale parameter in Gaussian kernel
K=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
KI=conv2(Img,K,'same');
KONE=conv2(ones(nrow, ncol),K,'same');
% start level set evolution
time = cputime;
penalizeType=0;
for n=1:iterNum
    method='ChunmingLi2005'; % this parameter specifies the method published in CVPR07 (developed in 2006).
    DiracFunction='global'; % this parameter specifies the Dirac function that is defined by the
                            % equation (12) in the CVPR07 paper. Do NOT change this parameter for this version. 
    Iter=1;
    % LSE: level set evolution.  
    u=LSE(u,Img,K,KI,KONE, nu,timestep,mu,lambda1,lambda2,epsilon,Iter,1,DiracFunction,method,2004);
    % The input '2004' represents the internal energy term used in my CVPR'05 and 07 papers (developed in 2004).
    % Do NOT change this parameter for this version.
    if mod(n,10)==0
        pause(0.001);
        imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
        contour(u,[0 0],'r');
        iterNum=[num2str(n*Iter), ' iterations'];
        title(iterNum);
        hold off;
    end
end
totaltime = cputime - time
imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
contour(u,[0 0],'r');
iterNum=[num2str(n), ' iterations'];
title(iterNum);


