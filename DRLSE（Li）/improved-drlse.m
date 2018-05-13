clear all;
close all;
%对李纯明的DRLSE进行了改进 利用各向异性扩散提高了分割弱边缘的能力 改进的高斯滤波代替惩罚项 加快演化速度 节省了时间
Img = imread('vessel.bmp'); % real miscroscope image of cells
Img=double(Img(:,:,1));
%% parameter setting
timestep=5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
%iter_inner=11;
iter_outer=100;
lambda=1.2; % coefficient of the weighted length term L(phi)
alfa=.8;  % coefficient of the weighted area term A(phi)
epsilon=1.5; % papramater that specifies the width of the DiracDelta function

sigma=1.5;     % scale parameter in Gaussian kernel
G=fspecial('gaussian',5,sigma);
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.

K1=fspecial('gaussian',9,sigma);  %%%%

% initialize LSF as binary step function
c0=2;
initialLSF=c0*ones(size(Img));
% generate the initial region R0 as a rectangle
initialLSF(5:55, 5:55)=-c0;  
%initialLSF(20:35, 35:55)=-c0; 
phi=initialLSF;

figure(1);
mesh(-phi);   % for a better view, the LSF is displayed upside down
hold on;  contour(phi, [0,0], 'r','LineWidth',2);
title('Initial level set function');
view([-80 35]);

figure(2);
imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
title('Initial zero level contour');
pause(0.5);



% start level set evolution
for n=1:iter_outer
   % phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
    [phi_x,phi_y]=gradient(phi);
    
    s=sqrt(phi_x.^2 + phi_y.^2);

    smallNumber=1e-10;  
    Nx=phi_x./(s+smallNumber); % add a small positive number to avoid division by zero
    Ny=phi_y./(s+smallNumber);
    %curvature=div(Nx,Ny);
    %f = div(nx,ny)
[nxx,junk]=gradient(Nx);  
[junk,nyy]=gradient(Ny);
curvature=nxx+nyy;
     %diracPhi=Dirac(phi,epsilon);
     f=(1/2/epsilon)*(1+cos(pi*phi/epsilon));
     b = (phi<=sigma) & (phi>=-sigma);
       diracPhi = f.*b;

     
     
    areaTerm=diracPhi.*g; % balloon/pressure force
    edgeTerm= diracPhi.*g.*curvature;%diracPhi.*(vx.*Nx+vy.*Ny) +
    phi=phi + timestep*(  lambda*edgeTerm + alfa*areaTerm);
    
    
    if mod(n,2)==0
        figure(2);
        imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
    end
    
    phi = (phi >= 0) - ( phi< 0);
    phi = conv2(phi, K1, 'same');
end

% refine the zero level contour by further level set evolution with alfa=0
%alfa=0;
%iter_refine = 10;
%phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);

finalLSF=phi;
figure(2);
imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
hold on;  contour(phi, [0,0], 'r');
str=['Final zero level contour, ', num2str(iter_outer), ' iterations'];
title(str);

pause(1);
figure;
mesh(-finalLSF); % for a better view, the LSF is displayed upside down
hold on;  contour(phi, [0,0], 'r','LineWidth',2);
str=['Final level set function, ', num2str(iter_outer), ' iterations'];
title(str);
axis on;

