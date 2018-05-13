function [u,f1,f2]= ACM(u,Img,Ksigma,KI,KI2,KONE,nu,timestep,mu,epsilon,lambda1,lambda2,CV,RSF,LRCV,LIF,LGD,isExchange)

u=NeumannBoundCond(u);
K=curvature_central(u);
Hu=0.5*(1+(2/pi)*atan(u./epsilon));
DrcU=(epsilon/pi)./(epsilon^2.+u.^2);

KIH= imfilter((Hu.*Img),Ksigma,'replicate');
KH= imfilter(Hu,Ksigma,'replicate');
f1=KIH./KH;
f2=(KI-KIH)./(KONE-KH);

[f1,f2]=exchange(f1,f2,isExchange); % important!

if CV~=0
    c= Hu.*Img;
    C1 = sum(c(:))/sum(Hu(:));
    c1=(1-Hu).*Img;
    c2=1-Hu;
    C2 = sum(c1(:))/sum(c2(:));
    CVterm=CV*(DrcU.*(-lambda1*(Img-C1).^2+lambda2*(Img-C2).^2));
else
    CVterm=0;
end

if LRCV~=0
    LRCVterm=LRCV*DrcU.*(-lambda1*(Img-f1).^2+lambda2*(Img-f2).^2);
else
    LRCVterm=0;
end

if LIF~=0
    LIFterm=DrcU.*((Img - f1.*Hu - f2.*(1 - Hu)).*(f1 - f2));
else
    LIFterm=0;
end

if RSF~=0
    s1=lambda1.*f1.^2-lambda2.*f2.^2;
    s2=lambda1.*f1-lambda2.*f2;
    dataForce=(lambda1-lambda2)*KONE.*Img.*Img+imfilter(s1,Ksigma,'replicate')-2.*Img.*imfilter(s2,Ksigma,'replicate');
    RSFterm=-RSF*DrcU.*dataForce;
else
    RSFterm=0;
end

if LGD~=0
    KI2H = imfilter(Img.^2.*Hu,Ksigma,'replicate');
    sigma1 = (f1.^2.*KH - 2.*f1.*KIH + KI2H)./(KH);
    sigma2 = (f2.^2.*KONE - f2.^2.*KH - 2.*f2.*KI + 2.*f2.*KIH + KI2 - KI2H)./(KONE-KH);
    localForce = (lambda1 - lambda2).*KONE.*log(sqrt(2*pi)) ...
        + imfilter(lambda1.*log(sqrt(sigma1)) - lambda2.*log(sqrt(sigma2)) ...
        +lambda1.*f1.^2./(2.*sigma1) - lambda2.*f2.^2./(2.*sigma2) ,Ksigma,'replicate')...
        + Img.*imfilter(lambda2.*f2./sigma2 - lambda1.*f1./sigma1,Ksigma,'replicate')...
        + Img.^2.*imfilter(lambda1.*1./(2.*sigma1) - lambda2.*1./(2.*sigma2),Ksigma,'replicate');
    LGDterm = -LGD*DrcU.*localForce;
else
    LGDterm=0;
end

PenaltyTerm=mu*(4*del2(u)-K);

LengthTerm=nu.*DrcU.*K;

u=u+timestep*(LengthTerm+PenaltyTerm+RSFterm+LRCVterm+LIFterm+LGDterm+CVterm);


function g = NeumannBoundCond(f)
% Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);

function k = curvature_central(u)
% compute curvature
[ux,uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);
Nx = ux./normDu;
Ny = uy./normDu;
[nxx,~] = gradient(Nx);
[~,nyy] = gradient(Ny);
k = nxx+nyy;

function [f1,f2]=exchange(f1,f2,isExchange)
%exchange f1 and f2
if isExchange==0
    return;
end
if isExchange==1
    f1=min(f1,f2);
    f2=max(f1,f2);
end
if isExchange==-1
    f1=max(f1,f2);
    f2=min(f1,f2);
end
