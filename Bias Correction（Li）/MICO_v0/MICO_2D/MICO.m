function [M_out, b_out, C_out]=  MICO(Img,q,W,M,C,b,Bas,GGT,ImgG, Iter, iterCM)
% This code implements the MICO algorithm for joint segmentation and  bias field estimation 
% proposed in the following paper:      
%      C. Li, J.C. Gore, and C. Davatzikos, 
%      "Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation", Magnetic Resonance
%      Imaging , vol. 32 (7), pp. 913-923, 2014
%
% All rights researved by Chunming Li
% E-mail: li_chunming@hotmail.com
% URL: http://imagecomputing.org/~cmli/
% Copyright (c) by Chunming Li
% Author: Chunming Li

for n = 1:Iter

    C = updateC(Img, W, b, M);
    for k=1:iterCM
        N_class=size(M,3);
        e=zeros(size(M));
        for kk=1:N_class
            D(:,:,kk) = (Img-C(kk)*b).^2;
        end
        M = updateM(D,q);        
    end
end

b_out = updateB(Img, q, C, M, Bas,GGT,ImgG);
M_out=M;
C_out=C;


%%%%%%%%%%%%%%%%%%
%%%  update b  %%%
function b =updateB(Img, q, C, M, Bas,GGT,ImgG)

PC2 = zeros(size(Img));
PC=PC2;

N_class=size(M,3);
for kk=1:N_class
    PC2=PC2+C(kk)^2*M(:,:,kk).^q;
    PC=PC+C(kk)*M(:,:,kk).^q;
end

N_bas=size(Bas,3);
V=zeros(N_bas,1);
for ii=1:N_bas
    ImgG_PC=ImgG{ii}.*PC;    % Mask in ImgG
    V(ii)=sum(ImgG_PC(:));   % inner product
    for jj=ii:N_bas
        B = GGT{ii,jj}.*PC2; % Mask in GGT
        A(ii,jj)=sum(B(:));  % inner product
        A(jj,ii)=A(ii,jj);
    end
end
clear PC1;
clear PC2;
clear B;
clear ImgG_PC;
w=inv(A)*V;

b=zeros(size(Img));
for kk=1:N_bas
    b=b+w(kk)*Bas(:,:,kk);
end

%%%%%%%%%%%%%%%%%%
%%%  update C  %%%
function C_new  =updateC(Img, W,b, M)
N_class=size(M,3);
for nn=1:N_class
    N=b.*Img.*M(:,:,nn);
    D=(b.^2) .*M(:,:,nn);
    sN = sum(N(:).*W(:));    % inner product
    sD = sum(D(:).*W(:));   % inner product
    C_new(nn)=sN/(sD+(sD==0));
end

clear N;
clear D;

%%%%%%%%%%%%%%%%%%
%%%  update M  %%%
function M = updateM(e, q)

N_class=size(e,3);

if q >1
    epsilon=0.000000000001;
    e=e+epsilon;  %% avoid division by zero
    p = 1/(q-1);
    f = 1./(e.^p);
    f_sum = sum(f,3);
    for kk=1:N_class
        M(:,:,kk) = f(:,:,kk)./f_sum;
    end
elseif q==1
    [e_min,N_min] = min(e,[], 3);  
    for kk=1:N_class
        M(:,:,kk) = (N_min == kk);
    end
else
    error('MICO: wrong fuzzifizer');
end


