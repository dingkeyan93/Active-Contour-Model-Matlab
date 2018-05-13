
clear all;close all;
im=imread('noisyImg.bmp');
im=double(im(:,:,1));

figure;imagesc(im, [0, 255]);colormap(gray);hold on;
[row,column]=size(im);
center=[row/2,column/2];
radius=20;
% initialize the phi function初始化水平集函数phi，其上各个网格点的值为该点到初始轮廓的距离
isinside=1;%
%phi = initphi( size( im ), center, radius, isinside );

BW = roipoly;   % get a region R inside a polygon, BW is a binary image with 1 and 0 inside or outside the polygon;
c0=4; % the constant value used to define binary level set function;
initialLSF= c0*2*(0.5-BW); % initial level set function: -c0 inside R, c0 outside R;
phi=initialLSF;
title('initial contour');
contour(phi,[0,0],'r');
% calculate the stopping function计算边缘停止函数
g = edgestop( im, 5, 2 );

% constants
delta_t  = 0.001;%时间步长
iterNum = 20;%边界轮廓迭代次数
bandIterNum = 10;%窄带的迭代次数
BAND_WIDTH = 2;%窄带宽度
iterNum_reinit=5;%iteration number of re-initialization

for n=1:iterNum

  % determine the front and narrow band points
  [ x, y ] = find( isfront( phi ) );%(仅)返回非零元素的行、列索引值,相当于一个二维数组
  front = [ x, y ];%边界点的索引值,
  %便于理解可以以如下形式考虑front[i][0]=x, front[i][1]=y.
  [ x, y ] = find( isband( phi, front, BAND_WIDTH ) );
  band = [ x, y ];%窄带内点的索引值

  phi=narrowbandEvolution(phi,g,front,band,delta_t,bandIterNum);
  %phi=reinit_SD_ENO2(phi, 1, 1, .5, iterNum_reinit); 
  pause(0.05);
  if mod(n,2)==0     
      imagesc(im,[0,255]);colormap(gray);hold on;
      [c,h] = contour(phi,[0 0],'r'); 
      iternum=[num2str(n), ' iterations'];        
      title(iternum);
      hold off;
  end;
  
end

imagesc(im,[0,255]);
colormap(gray);hold on;
[c,h]=contour(phi,[0,0],'g');
totalIterNum=['total iteration number is:',num2str(n)];
title(totalIterNum);      