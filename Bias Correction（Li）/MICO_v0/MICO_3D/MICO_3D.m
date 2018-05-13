% This function implements the MICO algorithm for joint segmentation and  bias field estimation 
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

function [M, b, C]=  MICO_3D(Img,W,M,C,b,Iter, iterCM,q)
for n = 1:Iter
    for kk=1:iterCM
        C = updateC(Img, W, q, b, M);   % Note: the order of updating C,M,and B can be changed.
        M = updateM(Img, W, C,b,q);
    end    
    b = updateB(Img, C, M,q);
end


%%%%%%%%%%%%%%%%%%
%%%  update b  %%%
function b =updateB(Img, C, M,q)

PC2 = zeros(size(Img));
PC=PC2;

N_class=size(M,4);
for kk=1:N_class
    Mq=M(:,:,:,kk).^q;
    PC2=PC2+C(kk)^2*Mq;
    PC=PC+C(kk)*Mq;
end

N_bas=20;

D=zeros(N_bas,1);
for i=1:N_bas
    filename = ['basis_',num2str(i),'.mat'];
    load(filename);
    basis_i = basis;
    A=Img.*basis_i.*PC;
    D(i)=sum(A(:));  
    for j=i:N_bas
        filename = ['basis_',num2str(j),'.mat'];
        load(filename);
        basis_j = basis;
        B = basis_i.*basis_j.*PC2;
        G(i,j)=sum(B(:));  
        G(j,i)=G(i,j);
    end
end
coeff=inv(G)*D;

clear A;
clear B;
clear PC;clear PC2;


b=zeros(size(Img));
for kk=1:N_bas
    filename = ['basis_',num2str(kk),'.mat'];
    load(filename);
    b=b+coeff(kk)*basis;
end

%%%%%%%%%%%%%%%%%%
%%%  update C  %%%
function C  =updateC(Img, W, q, b, M)
N_class=size(M,4);
for nn=1:N_class
    Mq=M(:,:,:,nn).^q;
    N=b.*Img.*Mq.*W;
    D=(b.^2 ).*Mq.*W;
    sN = sum(N(:));   
    sD = sum(D(:)); 
    C(nn)=sN/(sD+(sD==0));
end


%%%%%%%%%%%%%%%%%%
%%%  update M  %%%
function M_out = updateM(Img, W, C,b,q)

N_class=length(C);

if q > 1
    epsilon=0.000000000001;

    p = -1/(q-1);    
    f_sum = zeros(size(Img));
    for kk=1:N_class
        f_sum = f_sum + (((Img-C(kk)*b).^2) + epsilon).^p;
    end
    
    for kk=1:N_class
        M_out(:,:,:,kk) = W.*(((Img-C(kk)*b).^2) + epsilon).^p./f_sum;
    end
    clear f_sum;
    
elseif q==1
    for kk=1:N_class
        e(:,:,:,kk) = (Img-C(kk)*b).^2 ;
    end
    [e_min,N_min] = min(e,[], 4);   
    clear e;
    for kk=1:N_class
        M_out(:,:,:,kk) = (N_min == kk);
        M_out(:,:,:,kk) = M_out(:,:,:,kk).*W;
    end
    clear e N_min;
else
    error('MICO_3D: wrong fuzzifizer');
end

clear Img b;

