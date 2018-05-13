function u = levelset( im, center, radius, d_it, re_init, m_name )
%    Segment the image im using level set method, given 
%     an initial circle with the arguments: center and
%    radius. Then the current curve is displayed every d_it
%    iterations and written to disk s with name  m_name*.jpg


% set constants
delta_t  = 0.001;
ITERATIONS = 200000;
 
% initialize the u function
u = init_u( size( im ), center, radius, 1 );

% calculate the stopping function
g = stopfunction( im,3, 2 );

 [grad_g_x,grad_g_y]=gradient(g);


iterations = 0;

[I,J]=size(im);% image size

gauss_filter = fspecial( 'gaussian', 3, 2 );

while( iterations < ITERATIONS )
  if( mod( iterations, d_it ) == 0 )
    % display the segmented image every 'd_it' iterations
    % display current curve
    fprintf( 1, '%d\n', iterations );
    disp( 'Displaying segmented image' );
    segim = createimage( im, u );
    imshow( segim ); drawnow;
    
    filename = strcat( m_name, sprintf( '%06d', ( iterations / d_it ) + 1 ), '.jpg' );
    imwrite( segim, filename );
  end;
  % reinitialization of the embedding function u
  if( mod( iterations, re_init ) == 0 )
      u = re_init_u( u );
  end;
  % calculate the speed and update u
  speed = calcspeed( u, g, grad_g_x, grad_g_y ,iterations );

  u = u + (delta_t .* speed);
   
  if sum(sum(abs(delta_t .* speed)))<0.1
      break;
  end
  
  iterations = iterations + 1;    
end
% iterations