function seg = local_AC_UM(Img,mask_init,rad,alpha,num_it,epsilon)
% This function aims to implement Shawn Lankton's local active contour. And
% the energy model is the UM model as defined in eq.(11)-(14).
% The local variables are calculated by filtering operation instead of
% iterating inspired by Chunming Li's IEEE TIP 2008 paper
%
% One small change is that I used a square window instead of disk for
% localization
%
% Input: 
% 1. Img: image needs to be segmented
% 2. mask_init: intialization represented by binary image
% 3. rad: the side length of the square window
% 4. alpha: the coeficicent to balance the image fidality term and the
% curvature regularization term
% 5. num_it: maximum number of iterations
% 6. epsilon: epsilon used for delta and heaviside function
% Created by Jincheng Pang, Tufts University @11/09/2012

phi0 = mask2phi(mask_init);
phi = phi0;

B0 = ones(2*rad+1,2*rad+1);  
% B0 = fspecial('disk',rad);

KI=conv2(Img,B0,'same');  
KONE=conv2(ones(size(Img)),B0,'same'); 


for ii = 1:num_it
mask = Heaviside2(phi,epsilon);

I=Img.*mask;
temp1=conv2(mask,B0,'same');                             
temp2=conv2(I,B0,'same');                             
c1=temp2./(temp1);    % local mean value inside                                    
c2=(KI-temp2)./(KONE-temp1); % local mean value outside 

%%%%%% 
s1=(c1.^2-c2.^2).*Dirac2(phi,epsilon);                   % Compute the 1st term in eq. (14) in Lankton's IEEE TIP 08
s2=(c1-c2).*Dirac2(phi,epsilon);
dataForce=conv2(s1,B0,'same')-2.*Img.*conv2(s2,B0,'same'); %%% During the implementation, Img should be separated out of the filtering operation!!!
% dataForce=conv2(s1,B0,'same')-2.*conv2(s2.*Img,B0,'same'); %%% During the implementation, Img should be separated out of the filtering operation!!!
dataForce = dataForce/max(abs(dataForce(:)));
%%%%%%
% curvature = get_curvature1(phi);
curvature = curvature_central(phi); 
dphi = Dirac2(phi,epsilon).*(-dataForce + alpha*curvature);

dt = .48/(max(abs(dphi(:)))+eps);

%-- evolve the curve
phi = phi + dt.*dphi;

%-- Keep SDF smooth
 phi = sussman(phi, .5);

    if(mod(ii,10) == 0) 
      showCurveAndPhi(Img,phi,ii);  
    end
end

seg = (phi>=0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Auxiliary functions %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
  % forward/backward differences
  a = D - shiftR(D); % backward
  b = shiftL(D) - D; % forward
  c = D - shiftD(D); % backward
  d = shiftU(D) - D; % forward
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
  
%-- whole matrix derivatives
function shift = shiftD(M)
  shift = shiftR(M')';

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
  shift = shiftL(M')';
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);   

  function phi = mask2phi(init_a)
% Modified by Jincheng Pang to reverse the phi function;
%
%-- converts a mask to a SDF
  phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
  phi = -double(phi); % modified by Jincheng Pang 04/20/2012
  
  
  function showCurveAndPhi(I, phi, i)
  imshow(I,'initialmagnification',200,'displayrange',[ ]); 
% %   imagesc(xx,yy,I);axis square;axis xy
  hold on;  contour(phi, [0 0], 'y','LineWidth',2);
%   contour(phi, [0 0], 'k','LineWidth',4);
  hold off; title([num2str(i) ' Iterations']); drawnow;
  
  function f = Dirac2(x, sigma)
% % f=(1/2/sigma)*(1+cos(pi*x/sigma));
% % b = (x<=sigma) & (x>=-sigma);
% % f = f.*b;
f = (sigma/pi)./(sigma^2+x.^2);

  function f = Heaviside2(x, epsilon) % Use Heaviside_{2,epsilon} as denoted in Chan-Vese's TIP Paper.

     f = 0.5*(1+2/pi*atan(x./epsilon));

function k = curvature_central(u)                       
% compute curvature
[ux,uy] = gradient(u);                                  
normDu = sqrt(ux.^2+uy.^2+1e-10);                       % the norm of the gradient plus a small possitive number 
                                                        % to avoid division by zero in the following computation.
Nx = ux./normDu;                                       
Ny = uy./normDu;
[nxx,junk] = gradient(Nx);                              
[junk,nyy] = gradient(Ny);                              
k = nxx+nyy;                                            % compute divergence