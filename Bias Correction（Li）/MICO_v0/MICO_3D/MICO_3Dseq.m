% This Matlab function reads images in .nii format and calls the MICO
% function to perform simultaneous segmentation and bias field correction in Chunming Li et al's paper:
%    "Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation",
%     Magnetic Resonance Imaging, vol. 32 (7), pp. 913-923, 2014
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://imagecomputing.org/~cmli/
function MICO_3Dseq(str_vector, N_region, q, th_bg, iterNum_outer, Iter_b, iterCM, tissueLabel)

N_scan = length(str_vector);
for nn = 1:N_scan
    str=str_vector{nn};
    data = load_untouch_nii(str);
    Img=data.img;
    Img = double(Img);      
    save temp_info.mat data;
    clear data;    
    
    [x1, x2, y1, y2, z1, z2] = cropImg(Img, th_bg);   %% crop image
    [DimX1, DimY1, DimZ1]=size(Img);
    x1 = max(1,x1-2); x2 = min(DimX1,x2+2);
    y1 = max(1,y1-2); y2 = min(DimY1,y2+2);
    z1 = max(1,z1-2); z2 = min(DimZ1,z2+2);
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Img3D=Img(x1:x2,y1:y2,z1:z2);
    clear Img;
    [DimX, DimY, DimZ] = size(Img3D);    
    ROI = (Img3D>th_bg);
    saveBasisOrder3_3D(ROI);    
    Img3D = Img3D.*ROI;
    %%%%%%%%%%%%%%%%%%%%%% Initialization
    
    A=max(Img3D(:));
    C= linspace(0.1,0.9,N_region)*A;
    b=ones(size(Img3D));
    
    
    M=rand(DimX, DimY, DimZ,N_region);
    a=sum(M,4);
    for k = 1 : N_region
        M(:,:,:,k)=M(:,:,:,k)./a;
    end
    clear a;
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    totaltime = 0;
    M_old = M; chg=10000;
    save M_old M_old
    clear M_old;
   
    C_old =C;
    for n_iter= 1:iterNum_outer 
        n_iter
        tic
        [M, b, C]=   MICO_3D(Img3D,ROI,M,C,b,Iter_b,iterCM,q);
      
        totaltime = totaltime + toc
     
        [M, C]=sortMemC(M, C);
        PC2d=zeros(size(Img3D(:,:,1)));
        PC3d=zeros(DimX1, DimY1, DimZ1);
        N_slc=90;
        for k=1:N_region
            PC2d = PC2d +  tissueLabel(k)*M(:,:,N_slc,k);
        end
        pause(0.1);
        figure(1);
        subplot(1,2,1);
        imagesc(Img3D(:,:,N_slc));colormap(gray);
        subplot(1,2,2);
        imagesc(PC2d);colormap(gray);        
        C_new = C/norm(C);

        chg= max(abs(C_new(:)-C_old(:)));
        C_old = C_new;
        if chg<0.0001    % check convergence 
            break
        end         
    end   
    clear Img3D;
  
    %[U, C]=sortMemC(M, C);  
    U=maxMembership(M);
    clear M;
    Membership = zeros(DimX1, DimY1, DimZ1,N_region);    
    Membership(x1:x2,y1:y2,z1:z2,:)=U;  
    for k=1:N_region
        PC3d(:,:,:)=PC3d(:,:,:)+tissueLabel(k)*Membership(:,:,:,k);
    end
    pc3d=make_nii(PC3d);
    save_nii(pc3d,[str,'_seg.nii']);
    Bias = zeros(DimX1, DimY1, DimZ1);
    Bias(x1:x2,y1:y2,z1:z2)=b; clear b;
    bias=make_nii(Bias);
    save_nii(bias,[str,'_bc.nii']);

    
    clear Membership Bias U;
    
    for basis_index=1:20  % delete basis files
        filename = ['basis_',num2str(basis_index),'.mat'];
        delete(filename);
    end  
    delete M_old.mat temp_info.mat;
    
end

% sort the constants c1, c2, ..., and change the order of the membership
% functions accordingly.
function [M_out, C_out]=sortMemC(M, C)

[C_out IDX]=sort(C);

for k = 1 : length(C)
    M_out(:,:,:,k) = M(:,:,:,IDX(k));
end

% This Matlab function binarizes a fuzzy membership function 
function M_out = maxMembership(M)

if size(M,4)==1
    
    N_class=size(M,3);
    ROI = (sum(M,3)>0);
    
    M_out = zeros(size(M));
    
    [e_min,N_min] = max(M,[], 3);   % do not consider 2 minimum in this version
    for kk=1:N_class
        M_out(:,:,kk) = ROI.*(N_min == kk);
    end
    
elseif size(M,4)>1
    N_class=size(M,4);    
    ROI = (sum(M,4)>0);
  
        
    [e_min,N_min] = max(M,[], 4);   % do not consider 2 minimum in this version
    for kk=1:N_class
        M_out(:,:,:,kk) = ROI.*(N_min == kk);
    end
else
    error('wrong dimension: maxMembership');
end





