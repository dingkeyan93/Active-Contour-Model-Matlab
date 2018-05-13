% This Matlab file calls two localized active contour methods proposed by 
% Shawn Lankton's IEEE TIP 2008 paper 'Localizing Region-Based Active
% Contours'.
% function local_AC_UM is the localized Chan-Vese's method see TIP 2001.
% function local_AC_MS is the localized Antony's Mean Separation method see ICCV 1999  
%
% The image used in this script are from Chunming Li (http://www.engr.uconn.edu/~cmli/)

clc;clear all;close all;
imgID = 2; % 1,2,3  % choose one of the five test images

Img = imread([num2str(imgID),'.bmp']);
Img = double(Img(:,:,1));
epsilon = 1;
switch imgID

    case 1
        num_it =1000;
        rad = 8;
        alpha = 0.3;% coefficient of the length term
        mask_init  = zeros(size(Img(:,:,1)));
        mask_init(15:78,32:95) = 1;
        seg = local_AC_MS(Img,mask_init,rad,alpha,num_it,epsilon);
    case 2
        num_it =800;
        rad = 9;
        alpha = 0.003;% coefficient of the length term
        mask_init = zeros(size(Img(:,:,1)));
        mask_init(53:77,56:70) = 1;
        seg = local_AC_UM(Img,mask_init,rad,alpha,num_it,epsilon);
    case 3
        num_it = 1500;
        rad = 5;
        alpha = 0.001;% coefficient of the length term
        mask_init  = zeros(size(Img(:,:,1)));
        mask_init(47:80,86:99) = 1;
        seg = local_AC_UM(Img,mask_init,rad,alpha,num_it,epsilon);
end






