function f = sdf2circle(nrow,ncol, ic,jc,r)
%   sdf2circle(nrow,ncol, ic,jc,r) computes the signed distance to a circle
%   input: 
%       nrow: number of rows
%       ncol: number of columns
%       (ic,jc): center of the circle
%       r: radius of the circle
%   output: 
%       f: signed distance to the circle
%  
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li


[X,Y] = meshgrid(1:ncol, 1:nrow);

f = sqrt((X-jc).^2+(Y-ic).^2)-r;

%f=sdf2circle(100,50,51,25,10);figure;imagesc(f)