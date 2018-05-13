function g = edgestop( im, n, sigma )%定义一个边缘停止函数
% EDGESTOP Create the edge-stopping function
%    G = EDGESTOP( IM, SIZE, SIGMA ) creates the edge-stopping
%    function of 'im'. The image is first smoothed with a gaussian
%    kernal of size 'n' and parameter 'sigma'

% first convolve the image with a gaussian filter to blur it
gauss_filter = fspecial( 'gaussian', n, sigma );
filtered_im = imfilter( im, gauss_filter );
%函数imfilter:Multidimensional image filtering

% next create the gradient image by summing the horizontal
% and vertical gradient images
[ x_grad, y_grad ] = gradient( filtered_im );
grad_im = sqrt( ( x_grad.^2 ) + ( y_grad.^2 ) );

% scale the gradient image
max_grad = max( max( grad_im ) );
%括号内的max函数返回一个行向量，各个向量的值有各个列中的最大值组成。第二个max函数找到返回行向量中的最大值
min_grad = min( min( grad_im ) );%找到剃度平滑后图像中像素的最大和最小值
grad_im =  10 .* ( grad_im - min_grad ) ./ ( max_grad - min_grad );%scale the gradient image

% Create the edge-stopping function, can also use:
%g = exp( -abs( grad_im ) );
g = 1 ./ ( 1 + abs( grad_im ) );
