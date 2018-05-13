function speed = calcspeed( u, g, grad_g_x,grad_g_y,iterations)
% CALCSPEED Calculate the speed for all points

[I,J]=size(u);% image size

% fill speed matrix with zeros at first
speed = zeros( I , J );

% calculate the central differences
[Dx_central Dy_central] = gradient(u);

% calculate the forward differences
Dx_minus = u - circshift(u,[0,1]);
Dy_minus = u - circshift(u,[1,0]);

% calculate the backward differences
Dx_plus = circshift(u,[0,-1])-u;
Dy_plus = circshift(u,[-1,0])-u;

% calculate the curvature
u_x   = Dx_central;
u_y   = Dy_central;
u_xx  = Dx_plus - Dx_minus; %circshift(u,[0,-1]) + circshift(u,[0,1]) - 2 * u;
u_yy  = Dy_plus - Dy_minus; %circshift(u,[-1,0]) + circshift(u,[1,0]) - 2 * u;
u_xy1 = ( circshift(u,[-1,-1]) + circshift(u,[1,1]) ) / 4;
u_xy2 = ( circshift(u,[1,-1]) + circshift(u,[-1,1]) ) / 4;
u_xy  = u_xy1 - u_xy2;
K_top = u_xx.* ( u_y.^2 ) - ( 2 * u_x.*u_y.* u_xy ) + u_yy.* ( u_x.^2 );

K_bot = u_x.^2  + u_y.^2 ;
K_bot = max(K_bot,0.00001); 
F_remainder = K_top ./ K_bot ;

% estimate the final speed term
grad_g_u = grad_g_x .* Dx_central + grad_g_y .* Dy_central;                      
speed =g .*( F_remainder + 3 ) + grad_g_u; 