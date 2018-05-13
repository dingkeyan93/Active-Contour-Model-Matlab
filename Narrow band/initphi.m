function phi = initphi( imsize, center, radius, isinside )
% INITPHI Initialize the phi funcion
%    INITPHI( center, radius ) initializes the zero level level-set
%    of the phi function by the circle with the given center (x, y)
%    and radius. The size of phi is given by 'imsize'.
%用一个圆来初始化水平集函数phi的零水平集

% grab the requested size of phi
m = imsize( 1 ); n = imsize( 2 );%imsize=size(phi),imsize相当于一个数组,imsize(1)返回图像的行数。

% create a zero matrix
phi = zeros( imsize );

% go over each pixel
for i = 1 : m;
  for j = 1 : n;

    % get the sum of squares distance of the pixel from the center
    % of the circle　
    distance = sqrt( sum( ( center - [ i, j ] ).^2 ) );
    %distance为矩阵，其值为点(i,j)到中心点center的距离平方

    % set phi to be the signed distance from the pixel to the
    % circle's boundary, where the distance is positive for pixels
    % outside the boundary and negative for pixels inside the
    % boundary
    phi( i, j ) = distance - radius;

    if( isinside == 0 )
      phi( i, j ) = -phi( i, j );
    end

  end
end
