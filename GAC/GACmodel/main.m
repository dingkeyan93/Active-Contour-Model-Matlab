% 共10个文件(其中有1个图象"3.bmp",9个.m文件)
% 运行main函数,将对 3.bmp 图象进行分割.

% to open image file, set initial circle, and run either traditional segmentation 
clear all;
close all;

display_it = 100;
re_init = 1000;

% get the image filename and path
%[ filename, pathname ] = uigetfile( '*', 'Select an image to segment' );
filename = 'noisyNonUniform.bmp';

% if user cancelled dialog box then exit
if( filename == 0 )
  return;
end;

% Scale the image by a specified factor
scalefactor = 1;

% open the image and display it
% im = readimage( strcat( pathname, filename ));
im = readimage( filename );
im_resized = imresize( im, scalefactor );
imshow( uint8( im_resized ) );

% the circle center
center = floor(size(im_resized)/2);

%  define the radius
radius = min(center)-8;

moviename = filename( 1 : strfind( filename, '.' ) - 1 );

% run the level set segmentation
u =levelset(im_resized,circshift(center,[ 0, 1]),radius,display_it,re_init,moviename);