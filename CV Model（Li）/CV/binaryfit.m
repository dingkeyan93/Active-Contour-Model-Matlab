function [C1,C2]= binaryfit(Img,H_phi) 
%   [C1,C2]= binaryfit(phi,U,epsilon) computes c1 c2 for optimal binary fitting 
%   input: 
%       Img: input image
%       phi: level set function
%       epsilon: parameter for computing smooth Heaviside and dirac function
%   output: 
%       C1: a constant to fit the image U in the region phi>0
%       C2: a constant to fit the image U in the region phi<0
%  
%   Author: Chunming Li, all right reserved
%   email: li_chunming@hotmail.com
%   URL:   http://www.engr.uconn.edu/~cmli/research/

a= H_phi.*Img;
numer_1=sum(a(:)); 
denom_1=sum(H_phi(:));
C1 = numer_1/denom_1;

b=(1-H_phi).*Img;
numer_2=sum(b(:));
c=1-H_phi;
denom_2=sum(c(:));
C2 = numer_2/denom_2;
