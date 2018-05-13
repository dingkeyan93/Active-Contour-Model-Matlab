% This Matlab file demomstrates the method for simultaneous segmentation and bias field correction 
% in Chunming Li et al's paper:
%    "Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation",
%     Magnetic Resonance Imaging, vol. 32 (7), pp. 913-923, 2014
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://imagecomputing.org/~cmli/
clear all;
close all;

iterNum_outer=15;  % outer iteration
iterCM=2;  % inner interation for C and M
iter_b=1;  % inner iteration for bias

A = 255;
q = 1.5;

tissueLabel=[1, 2, 3];

th_bg = 5;  %% threshold for removing background
N_region = 3; %% number of tissues, e.g. WM, GM, CSF

str_vector{1} = 'brainweb_byte_B3N2.nii';   % input a sequence of image file names
% str_vector{2} = 'brainweb_byte_B3N1.nii';


MICO_3Dseq(str_vector, N_region, q, th_bg, iterNum_outer, iter_b, iterCM, tissueLabel);
