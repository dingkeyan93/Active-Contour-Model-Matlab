% This Matlab file demomstrates a narrow band algorithm that implements the level set method in Li et al's paper
%    "Level Set Evolution Without Re-initialization: A New Variational Formulation"
%    in Proceedings of CVPR'05, vol. 1, pp. 430-436.
% Author: Chunming Li, all rights reserved.
% E-mail: li_chunming@hotmail.com
% URL:  http://www.engr.uconn.edu/~cmli/

clear all;
close all;
Img = imread('twoObjImg.bmp');
Img=double(Img(:,:,1));
sigma=1.5;    % scale parameter in Gaussian kernel for smoothing.
[nrow, ncol]=size(Img);
total_time = 0;

G=fspecial('gaussian',15,sigma);
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.
[vx,vy] = gradient(g);
tic;
total_time = total_time + toc;
epsilon=1.5; % the papramater in the definition of smoothed Dirac function

N_iter=150;
timestep=5;  % time step
mu=0.04;  % coefficient of the internal (penalizing) energy term P(\phi)
% Note: The product timestep*mu must be less than 0.25 for stable evolution
DiracSigma=1.5; % the papramater in the definition of smoothed Dirac function
lambda=5; % coefficient of the weighted length term Lg(\phi)
alf=3;   % coefficient of the weighted area term Ag(\phi);
% Note: Choose a positive(negative) alf if the initial contour is outside(inside) the object.


d = 1;  % specify the width of the narrow band. The narrow band is
% the union of (2d+1)x(2d+1) blocks centered at the zero crossing pixels.
% Note: don't change this parameter for this version (nband_v14).
epsilon = .000001;  % small positive number added to the denominator to avoid 0 in the denominator
rad = d;
vs=d+1;

figure;imagesc(Img, [0, 255]);colormap(gray);hold on;
text(6,6,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'r');
% Use mouse to obtain initial contour/region;
BW = roipoly;  % get a region R inside a polygon, BW is a binary image with 1 and 0 inside or outside the polygon;
c0=2;  % The constant value used to define binary level set function as initial LSF;
% Using larger value of c0 usually slows down the evolution.

initialLSF= c0*2*(0.5-BW); % initial level set function: -c0 inside R, c0 outside R;
u=initialLSF;
figure;imagesc(Img, [0, 255]);colormap(gray);hold on;
[c,h] = contour(u,[0 0],'b');
title('Initial contour');

%initial bandmap,
bandMap = initializeNB(u);

total_iter_time = 0;

ops{1}.name     = 'normalGradient';
ops{end+1}.name = 'Dirac';
ops{end+1}.name = 'curvature';
ops{end+1}.name = 'computeDataTerm';
ops{end+1}.name = 'updateLSF';
ops{end+1}.name = 'updateZeroCrossBand';
ops{end+1}.name = 'updateZCBand init';
ops{end+1}.name = 'total time';

total_stats = zeros(length(ops), 1);

tic;
t1=cputime;
[u, bandMap, ux, uy, Nx, Ny, diracU, K, ABD, f, E,Stats] = ...
    nband_v14(u, bandMap, DiracSigma, vs, lambda, alf, mu, timestep, rad, vx, vy, g, nrow, ncol, N_iter);
t2=cputime;
iterationTime=t2-t1
total_stats = total_stats + Stats;
total_iter_time = total_iter_time + toc;
total_time = total_time + total_iter_time;

%hold on;
figure;imagesc(Img, [0, 255]);colormap(gray);hold on;
[c,h] = contour(u,[0 0],'r');
iterNum=[num2str(N_iter), ' iterations'];
title(iterNum);


disp(' ');
printStats(ops, total_stats);
disp(' ');

% display statistics of CPU times
disp('  =============== time statistics (Matlab) ================');
disp(sprintf('%26s\t %.6g s', 'total iteration time', total_iter_time));
disp(sprintf('%26s\t %.6g s', 'total time (init/iter)', total_time));
disp(' ');




