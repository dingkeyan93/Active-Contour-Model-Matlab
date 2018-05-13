function [c,h]=plotLevelSet(u,zLevel, style)
%   plotLevelSet(u,zLevel, style) plot the level contour of function u at
%   the zLevel.
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li
% hold on;
[c,h] = contour(u,[zLevel zLevel],style); 
% hold off;