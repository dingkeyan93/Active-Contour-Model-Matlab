function [c,h]=plotLevelSet(u,zLevel, style)
%   plotLevelSet(u,zLevel, style) plot the level contour of function u at

[c,h] = contour(u,[zLevel zLevel],style); 
