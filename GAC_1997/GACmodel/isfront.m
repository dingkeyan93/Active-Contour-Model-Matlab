function front = isfront(u)
% ISFRONT Determine whether a pixel is a front point
%    ISFRONT(u ) return binary matrix whose value at each pixel
%    represents whether the corresponding pixel in u is a front
%    point or not.

% grab the size ofu
[ n, m ] = size(u );

% create an boolean matrix whose value at each pixel is 0 or 1
% depending on whether that pixel is a front point or not
front = zeros( size(u ) );

% A piecewise linear approximation to the front is contructed by
% checking each pixels neighbour. Do not check pixels on border.
for i = 2 : n - 1;
  for j = 2 : m - 1;

    % if there is a sign change then we have a front point
    maxVal = max( max(u( i:i+1, j:j+1 ) ) );
    minVal = min( min(u( i:i+1, j:j+1 ) ) );
    front( i, j ) = ( ( maxVal > 0 ) & ( minVal < 0 ) ) | u( i, j ) == 0;
 
  end
end