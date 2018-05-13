function front = isfront( phi )%返回一个二值矩阵，矩阵的每个像素点的值表示当前点是否为边界点
% ISFRONT Determine whether a pixel is a front point
%    ISFRONT( phi ) return binary matrix whose value at each pixel
%    represents whether the corresponding pixel in phi is a front
%    point or not.

% grab the size of phi
[ n, m ] = size( phi );

% create an boolean matrix whose value at each pixel is 0 or 1
% depending on whether that pixel is a front point or not
front = zeros( size( phi ) );

% A piecewise(分段) linear approximation to the front is contructed by
% checking each pixels neighbour. Do not check pixels on border.
for i = 2 : n - 1;
  for j = 2 : m - 1;

    % if there is a sign change then we have a front point
    %判断(i,j),(i+1,j),(i,j+1)和(i+1,j+1)四个点中是否包含有边界点
    %****注意：求矩阵元素的最大或最小值的形式为：max(max(matrix))或min(min(matrix))
    maxVal = max( max( phi( i:i+1, j:j+1 ) ) );%找出四个点中的最大值
    minVal = min( min( phi( i:i+1, j:j+1 ) ) );%找出四个点中的最小值
    front( i, j ) = ( ( maxVal > 0 ) & ( minVal < 0 ) ) | phi( i, j ) == 0;
    %当四个点中的最大和最小值异号或者当前点(i,j)值为零，则包含有边界点

  end
end