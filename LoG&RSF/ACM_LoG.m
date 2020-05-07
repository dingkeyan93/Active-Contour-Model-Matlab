function [u,f1,f2]= ACM_LoG(u,Img,Ksigma,KI,KONE,nu,timestep,mu,epsilon,lambda1,lambda2,RSF,LRCV,LIF,isExchange,log_)

u=NeumannBoundCond(u);
K=curvature_central(u);

Hu=0.5*(1+(2/pi)*atan(u./epsilon));
DrcU=(epsilon/pi)./(epsilon^2.+u.^2);

KIH= imfilter((Hu.*Img),Ksigma,'replicate');
KH= imfilter(Hu,Ksigma,'replicate');
f1=KIH./KH;
f2=(KI-KIH)./(KONE-KH);

[f1,f2]=exchange(f1,f2,isExchange);

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

LoGterm=log_.*DrcU; %

PenaltyTerm=mu*(4*del2(u)-K);

LengthTerm=nu.*DrcU.*K;

u = u + timestep*(LengthTerm+PenaltyTerm+RSFterm+LRCVterm+LIFterm+LoGterm); % update level set


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
