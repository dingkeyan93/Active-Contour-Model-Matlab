function K = CURVATURE(f,diff_scheme)
% CURVATURE computes curvature
% Author: Chunming Li, all rights reserved.
% Email: li_chunming@hotmail.com
% URL: http://www.engr.uconn.edu/~cmli/

epsilon=1e-10;

if strcmp(diff_scheme, 'fcb')  
    [fx,fy]=gradient(f);  % central difference
    fx_f = Dx_forward(f); % forward difference
    ax = fx_f./sqrt(fx_f.^2+ fy.^2+epsilon);
    axx = Dx_backward(ax); % backward difference
    fy_f = Dy_forward(f);
    ay = fy_f./sqrt(fx.^2 + fy_f.^2 + epsilon);
    ayy = Dy_backward(ay); 
    K = axx + ayy;

elseif strcmp(diff_scheme, 'fb')   % forward difference followed by a backward difference

    fx_f = Dx_forward(f);
    fy_f = Dy_forward(f);
    ax = fx_f./sqrt(fx_f.^2+ fy_f.^2+epsilon);
    ay = fy_f./sqrt(fx_f.^2 + fy_f.^2 + epsilon);
    axx = Dx_backward(ax);
    ayy = Dy_backward(ay);
    K = axx + ayy;
elseif strcmp(diff_scheme, 'bf')   % forward difference followed by a backward difference

    fx_f = Dx_backward(f);
    fy_f = Dy_backward(f);
    ax = fx_f./sqrt(fx_f.^2+ fy_f.^2+epsilon);
    ay = fy_f./sqrt(fx_f.^2 + fy_f.^2 + epsilon);
    axx = Dx_forward(ax);
    ayy = Dy_forward(ay);
    K = axx + ayy;
elseif strcmp(diff_scheme, 'cc')   % central difference followed by a central difference
    [fx, fy]= gradient(f); % central difference
    ax = fx./sqrt(fx.^2+ fy.^2+epsilon);
    ay = fy./sqrt(fx.^2 + fy.^2 + epsilon);
    [axx, axy] = gradient(ax); % central difference
    [ayx, ayy] = gradient(ay);
    K = axx + ayy;    
else
    disp('Wrong difference scheme: CURVATURE.m');
    return;    
end





