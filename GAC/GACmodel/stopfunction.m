function g = stopfunction( im, n, sigma )

%    This function calculates the edge-stopping
%    function of 'im'. The 'im' is firstly smoothed with a gaussian
%    kernal of size 'n' and parameter 'sigma', and then calculate its
%    gradeint vector

% first convolve the image with a gaussian filter 
gauss_filter = fspecial( 'gaussian', n, sigma );
filtered_im = imfilter( im, gauss_filter );
[ x_grad, y_grad ] = gradient( filtered_im );

grad_im = sqrt( ( x_grad.^2 ) + ( y_grad.^2 ) );

%  normolize the gradient image
max_grad = max( max( grad_im ) );
min_grad = min( min( grad_im ) );
grad_im =  100 .* ( grad_im - min_grad ) ./ ( max_grad - min_grad );

% Create the edge-stopping function:
g = exp( -abs( grad_im ) );
%g = 1 ./ (( 1 + abs( grad_im )).^2 );

