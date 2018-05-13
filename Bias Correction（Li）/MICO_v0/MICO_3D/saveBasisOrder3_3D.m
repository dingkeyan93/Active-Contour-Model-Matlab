function saveBasisOrder3_3D(ROI)

[DimX, DimY,DimZ]=size(ROI);
% x = zeros(Height,Wide,Z);
% y = zeros(Height,Wide,Z);
% z = zeros(Height,Wide,Z);

tempX = -1:2/(DimX-1):1;
tempY = -1:2/(DimY-1):1;
tempZ = -1:2/(DimZ-1):1;

[x y z]=meshgrid(tempY, tempX, tempZ);
basis=ones(size(x));
for basis_index = 1:20
    
    switch basis_index
        case 1
            basis = ROI;
        case 2
            basis = ROI.*z;
        case 3
            basis = ROI.*(3*z.*z - 1)/2;
        case 4
            basis = ROI.*(5*z.*z.*z - 3*z)/2;
        case 5
            basis = ROI.*y;
        case 6
            basis = ROI.*y.*z;
        case 7
            basis = ROI.*y.*(3*z.*z - 1)/2;
            
        case 8
            basis = ROI.*(3*y.*y - 1)/2;
        case 9
            basis = ROI.*z.*(3*y.*y - 1)/2;
            
        case 10
            basis = ROI.*(5*y.*y.*y - 3*y)/2;
            
        case 11
            basis = ROI.*x;
        case 12
            basis = ROI.*x.*z;
        case 13
            basis = ROI.*x.*(3*z.*z-1)/2;
            
        case 14
            basis = ROI.*x.*y;
        case 15
            basis = ROI.*x.*y.*z;
            
        case 16
            basis = ROI.*x.*(3*y.*y-1)/2;
            
        case 17
            basis = ROI.*(3*x.*x-1)/2;
        case 18
            basis = ROI.*z.*(3*x.*x-1)/2;
            
        case 19
            basis = ROI.*y.*(3*x.*x-1)/2;
            
        case 20
            basis = ROI.*(5*x.*x.*x-3*x)/2;
    end
  
    A=basis.^2;
    s = sum(A(:));
    basis=basis/s;
    clear A s;

    filename = ['basis_',num2str(basis_index),'.mat'];
    save(filename,'basis');
    
end

clear x y z basis;

