function u = init_u( imsize, center, radius, isinside )
%  init_u Initialize the u funcion
%    initializes the zero level level-set
%    of the u function by the circle with the given center (x, y)
%    and radius. The size of u is given by 'imsize'.

% grab the requested size of u
m = imsize( 1 ); n = imsize( 2 );

% create a zero matrix
u = zeros( imsize );

% go over each pixel
for i = 1 : m;
  for j = 1 : n;

    % get the sum of squares distance of the pixel from the center
    % of the circle
    distance = sqrt( sum( ( center - [ i, j ] ).^2 ) );

    % set u to be the signed distance from the pixel to the
    % circle's boundary, where the distance is positive for pixels
    % outside the boundary and negative for pixels inside the
    % boundary
    u( i, j ) = distance - radius;

    % if( isinside == 0 )
    %  u( i, j ) = -u( i, j );
    % end

  end
end
