function newim = createimage( im, u )
% CREATEIMAGE Return segmented image as 3 channel uint8 image
%    CREATEIMAGE( im, u ) Draws the segmented region in red
%    overtop of the original image im. The region is determined by
%    finding the front points from u

% copy over pixels to a temporary image and outline pixels of the
% segmented region
tempim = im;
% copy over pixels into red channel
front = isfront( u );
front = circshift(front,[1,0]);
tempim( find( front ) ) = 255;
newim( :, :, 1 ) = tempim;
tempim( find( front ) ) = 0;
newim( :, :, 2 ) = tempim;
tempim( find( front ) ) = 0;
newim( :, :, 3 ) = tempim;
% convert to uint8
newim = uint8( newim );
