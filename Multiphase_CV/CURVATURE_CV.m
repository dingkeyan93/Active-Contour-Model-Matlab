function K = CURVATURE_CV(f,diff_scheme)
% K = CURVATURE_CV(f,diff_scheme) computes curvature
% Copyright (c) 2001--2006 by Chunming Li
% Author: Chunming Li, 11/02/2003
% Revision by Chunming Li 8/18/2006

epsilon=1e-10;

if diff_scheme == 0  % the scheme in Chan and Vese's paper
    [fx,fy]=gradient(f);  % central difference
    fx_f = forward_gradient(f);
    ax = fx_f./sqrt(fx_f.^2+ fy.^2+epsilon);
    axx = backward_gradient(ax);    
    fy_f = forward_gradient(f);
    ay = fy_f./sqrt(fx.^2 + fy_f.^2 + epsilon);
    ayy = backward_gradient(ay);    
    K = axx + ayy;
    
elseif diff_scheme == 1   % forward difference followed by a backward difference
    
    fx_f = Dx_forward(f);
    fy_f = Dy_forward(f);
    ax = fx_f./sqrt(fx_f.^2+ fy_f.^2+epsilon);
    ay = fy_f./sqrt(fx_f.^2 + fy_f.^2 + epsilon);    
    axx = Dx_backward(ax);
    ayy = Dy_backward(ay);    
    K = axx + ayy;
elseif diff_scheme == 2   % central difference followed by a central difference    
    [fx, fy]= gradient(f); % central difference
    ax = fx./sqrt(fx.^2+ fy.^2+epsilon);
    ay = fy./sqrt(fx.^2 + fy.^2 + epsilon);    
    [axx, axy] = gradient(ax); % central difference
    [ayx, ayy] = gradient(ay);    
    K = axx + ayy;    
else
    disp('Wrong difference scheme: CURVATURE_CV.m');
    return;    
end





