function H = Heaviside(phi,epsilon) 
%   Heaviside(phi,epsilon)  compute the smooth Heaviside function
%  
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li
H = 0.5*(1+ (2/pi)*atan(phi./epsilon));