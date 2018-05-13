function [x1 x2 y1 y2 z1 z2] = cropImg(Img, th)

[Nx, Ny, Nz]=size(Img);

xx=[];
yy=[];
zz=[];


x1=Nx;
x2=1;
y1=Ny;
y2=1;
z1=Nz;
z2=1;

zz=[];

for z=1:Nz
    Img2D=Img(:,:,z);
    [xx,yy]=find(Img2D>th);
    if length(xx)>0
        zz=[zz,z];
        x1=min(x1,min(xx));
        x2=max(x2,max(xx));
        y1=min(y1, min(yy));
        y2=max(y2, max(yy));
    end



end

z1=min(zz);
z2=max(zz);





