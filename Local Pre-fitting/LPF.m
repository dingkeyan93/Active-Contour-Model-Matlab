function [f1,f2,s1,s2,Im] = LPF(Img,K)
[m2,n2]=size(K);
KK=K*m2*n2;
r=(m2-1)/2;
[m,n]=size(Img);
imgn=zeros(m+2*r,n+2*r);
imgn(r+1:m+r,r+1:n+r)=Img;
imgn(1:r,r+1:n+r)=Img(1:r,1:n);                 %扩展上边界
imgn(1:m+r,n+r+1:n+2*r+1)=imgn(1:m+r,n:n+r);    %扩展右边界
imgn(m+r+1:m+2*r+1,r+1:n+2*r+1)=imgn(m:m+r,r+1:n+2*r+1);    %扩展下边界
imgn(1:m+2*r+1,1:r)=imgn(1:m+2*r+1,r+1:2*r);       %扩展左边界
f1=zeros(m+2*r,n+2*r);
f2=zeros(m+2*r,n+2*r);
s1=zeros(m+2*r,n+2*r);
s2=zeros(m+2*r,n+2*r);
Im=zeros(m+2*r,n+2*r);
for i=r+1:m+r
    for j=r+1:n+r
        tmpM=imgn(i-r:i+r,j-r:j+r);
        Im(i,j)=mean(tmpM(tmpM~=0));
        ind1=find(tmpM>0&tmpM<=Im(i,j));
        ind2=find(tmpM>=Im(i,j));
        tmpM=tmpM.*KK;
        f1(i,j)=sum(tmpM(ind1))/(sum(KK(ind1))+eps);
        s1(i,j)=sum((f1(i,j)-tmpM(ind1)).^2);
        f2(i,j)=sum(tmpM(ind2))/(sum(KK(ind2))+eps);
        s2(i,j)=sum((f2(i,j)-tmpM(ind2)).^2);
    end
end
f1=f1(r+1:m+r,r+1:n+r);
f2=f2(r+1:m+r,r+1:n+r);
s1=s1(r+1:m+r,r+1:n+r);
s2=s2(r+1:m+r,r+1:n+r);
Im=Im(r+1:m+r,r+1:n+r);

