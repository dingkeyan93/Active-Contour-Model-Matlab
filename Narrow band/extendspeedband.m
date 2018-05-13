function speed = extendspeedband( speed, front, band )
% EXTENDSPEEDBAND Extends speed to non-front points扩展速度到窄带内的非边缘点
%    EXTENDSPEEDBAND( speed, front ) Extends speed function to
%    all pixels in a narrow band around 'front' given by the
%    indices in 'band'. For each band pixel, its speed is
%    the same as the speed of the closest front　pixel.
%    扩展速度函数到窄带内的所有点，每个窄带内点的速度和距离它最近的边缘上的点的速度相同

% grab the size of the speed matrix
[ m, n ] = size( speed );

% grab the total number of front and band points
n_front = size( front, 1 );%边缘点的数目
n_band  = size( band, 1 );%窄带内点的数目

% for every pixel in the band
for ii = 1 : n_band;

  % grab the coordinates of the band pixel获得窄带内点的坐标索引值
  i = band( ii, 1 ); j = band( ii, 2 );

  % find the front pixel it is closest to找到距离点（i,j）最近的边缘点
  closest_dist = inf;
  closest_point = front( 1, : );%front(k,:)返回第ｋ个点的x,y坐标值
  for kk = 1 : n_front;
    dist = sum( ( front( kk, : ) - [ i, j ] ).^2 );
    if( dist < closest_dist )
      closest_dist = dist;
      closest_point = front( kk, : );
    end;
  end;

  % and copy its speed
  speed( i, j ) = speed( closest_point( 1 ), closest_point( 2 ) );
  %closest_point(1)返回x坐标的索引值,closest_point(2)返回ｙ坐标索引值
end;
