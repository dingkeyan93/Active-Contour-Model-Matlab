function speed = extendspeed( speed, front )%扩展速度到全局内的非边界点（非窄带）
% EXTENDSPEED Extends speed to non-front points
%    EXTENDSPEED( speed, front ) Extends speed function to
%    all pixels in 'speed'. For each non-front pixel, its speed is
%    the same as the speed of the closest front pixel.

% grab the size of the speed matrix
[ m, n ] = size( speed );

% grab the total number of front points
n_front = size( front, 1 );

% for every pixel找出点(i,j)在边界线上的最近点P，并把P点的速度赋给点(i,j)
for i = 1 : m;
  for j = 1 : n;

    % find the front pixel it is closest to
    closest_dist = inf;
    closest_pt = front( 1, : );
    for k = 1 : n_front;
      dist = sum( ( front( k, : ) - [ i, j ] ).^2 );
      if( dist < closest_dist )
	closest_dist = dist;
	closest_pt = front( k, : );%把第ｋ个点，即最近点，赋值给closest_pt.closest_pt(1)表示点的ｘ坐标
      end;
    end;

    % and copy its speed
    speed( i, j ) = speed( closest_pt( 1 ), closest_pt( 2 ) );
  end;
end;
