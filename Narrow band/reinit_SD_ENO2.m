function [phi] = reinit_SD_ENO2(phi, dx, dy, alpha, iterations)
%该方法见文件夹中的PPT"Level Set方法及其在图像图形中的应用"
% function [phi] = reinit_SD(phi, dx, dy, alpha, accuracy, iterations)
%重新初始化phi为符号距离函数
% Reinitializes phi into a signed distance function while preserving
% the zero level set (the interface or the curve).
%
% dx and dy are the resolution of the grid at x and y dimensions.
% alpha is a constant for calculating the euler step (dt). Should
% be between 0 and 1. 0.5 is quite safe whereas 0.9 can be risky.
% iterations specifies the number of iterations before the function returns.
% accuracy is the order of accuracy of derivative calculation.
% Allowed values for accuracy are 'ENO1', 'ENO2', 'ENO3', 'WENO'. 
% These correspond to 1st, 2nd, 3rd and 5th order accurate schemes 
% for calculating the derivative of phi.
%alpha为重新初始化的时间步长
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/bwdist


init_normal = @init_normal_ENO2;
evolve_normal = @evolve_normal_ENO2;

S_phi_0 = phi./sqrt(phi.^2 + dx.^2);

Vn_ext = feval(init_normal, S_phi_0);
it=0;
t=0;
while(it < iterations)
	[delta_normal, H1_abs, H2_abs] = feval(evolve_normal, phi, dx, dy, Vn_ext);
%feval(' fun ',x)求由字符串' fun '给定的函数值，其输入参量是变量x。即feval(' fun ',x)等价于求fun(x)值。

	dt = get_dt_normal(alpha, dx, dy, H1_abs, H2_abs);
	phi = phi + dt*(S_phi_0 - delta_normal);
	it = it+1;
	t = t+dt;
end


function [dt] = get_dt_normal(alpha, dx, dy, H1_abs, H2_abs)%计算导数的法向
if alpha <= 0 | alpha >= 1 
    error('alpha needs to be between 0 and 1!');
end

maxs = max(max(H1_abs/dx + H2_abs/dy));
dt = alpha/(maxs+(maxs==0));

function [Vn_ext] = init_normal_ENO2(Vn)
Vn_ext = zeros(size(Vn)+4);
Vn_ext(3:end-2,3:end-2) = Vn;

function [der] = select_der_normal(Vn, der_minus, der_plus)

if size(der_minus) ~= size(der_plus) | size(der_plus) ~= size(Vn)
    error('plus, minus derivative vectors and normal force (Vn) need to be of equal length!');
end

der = zeros(size(der_plus));

for i=1:numel(Vn)
	Vn_der_m = Vn(i)*der_minus(i);
	Vn_der_p = Vn(i)*der_plus(i);
	if Vn_der_m <= 0 & Vn_der_p <= 0
		der(i) = der_plus(i);
	elseif Vn_der_m >= 0 & Vn_der_p >= 0
		der(i) = der_minus(i);
	elseif Vn_der_m <= 0 & Vn_der_p >= 0
		der(i) = 0;
	elseif Vn_der_m >= 0 & Vn_der_p <= 0
		if abs(Vn_der_p) >= abs(Vn_der_m)
			der(i) = der_plus(i);
		else
			der(i) = der_minus(i);
		end
	end
end


function [data_x] = der_ENO2_minus(data, dx)%dx、dy分别为x,y轴的像素分辨率
data_x = zeros(size(data));

% extrapolate the beginning and end points of data对起始点的像素进行插值
data(2) = 2*data(3)-data(4);
data(1) = 2*data(2)-data(3);
data(end-1) = 2*data(end-2)-data(end-3);
data(end) = 2*data(end-1)-data(end-2);

%Generate the divided difference tables(均差, 除法差分表)
%ignoring division by dx for efficiency
D1 = (data(2:end)-data(1:end-1));%向前差分
% ignoring division by dx since this will cancel out(抵偿)
D2 = (D1(2:end)-D1(1:end-1))/2;%二次插值
absD2 = abs(D2);%求绝对值

for i=1:(length(data)-4)
%length(data)返回矩阵data的行数或列数中的最大值. 
    k = i-1;

    Q1p = D1(k+2); %D1k_half;

    if absD2(k+1) <= absD2(k+2) %|D2k| <= |D2kp1|
        c = D2(k+1); %D2k;
    else
        c = D2(k+2); %D2kp1;
    end

    % ignoring multiplication by dx since this will also cancel out
    Q2p = c*(2*(i-k)-1);

    data_x(i+2) = Q1p+Q2p;
    data_x(i+2) = data_x(i+2)/dx;
end


function [data_x] = der_ENO2_plus(data, dx)
data_x = zeros(size(data));

% extrapolate the beginning and end points of data
data(2) = 2*data(3)-data(4);%对新增加的边缘行、列进行插值
data(1) = 2*data(2)-data(3);
data(end-1) = 2*data(end-2)-data(end-3);
data(end) = 2*data(end-1)-data(end-2);

%Generate the divided difference tables
%ignoring division by dx for efficiency
D1 = (data(2:end)-data(1:end-1));
% ignoring division by dx since this will cancel out
D2 = (D1(2:end)-D1(1:end-1))/2;
absD2 = abs(D2);

for i=1:(length(data)-4)
    k = i;

    Q1p = D1(k+2); %D1k_half;

    if absD2(k+1) <= absD2(k+2) %|D2k| <= |D2kp1|
        c = D2(k+1); %D2k;
    else
        c = D2(k+2); %D2kp1;
    end

    % ignoring multiplication by dx since this will also cancel out
    Q2p = c*(2*(i-k)-1);

    data_x(i+2) = Q1p+Q2p;
    data_x(i+2) = data_x(i+2)/dx;
end


function [delta, H1_abs, H2_abs] = evolve_normal_ENO2(phi, dx, dy, Vn)
delta = zeros(size(phi)+4);%二阶精度，所以行、列要分别加上４。三阶时则要加６
data_ext = zeros(size(phi)+4);
data_ext(3:end-2,3:end-2) = phi;

% Calculate the derivatives (both + and -)
phi_x_minus = zeros(size(phi)+4);
phi_x_plus = zeros(size(phi)+4);
phi_y_minus = zeros(size(phi)+4);
phi_y_plus = zeros(size(phi)+4);
phi_x = zeros(size(phi)+4);
phi_y = zeros(size(phi)+4);
% first scan the rows
for i=1:size(phi,1)%size(phi,1)返回矩阵phi的行数
	phi_x_minus(i+2,:) = der_ENO2_minus(data_ext(i+2,:), dx);	
	phi_x_plus(i+2,:) = der_ENO2_plus(data_ext(i+2,:), dx);	
	phi_x(i+2,:) = select_der_normal(Vn(i+2,:), phi_x_minus(i+2,:), phi_x_plus(i+2,:));
end

% then scan the columns
for j=1:size(phi,2)
	phi_y_minus(:,j+2) = der_ENO2_minus(data_ext(:,j+2), dy);	
	phi_y_plus(:,j+2) = der_ENO2_plus(data_ext(:,j+2), dy);	
	phi_y(:,j+2) = select_der_normal(Vn(:,j+2), phi_y_minus(:,j+2), phi_y_plus(:,j+2));
end

abs_grad_phi = sqrt(phi_x.^2 + phi_y.^2);

H1_abs = abs(Vn.*phi_x.^2 ./ (abs_grad_phi+dx*dx*(abs_grad_phi == 0)));
H2_abs = abs(Vn.*phi_y.^2 ./ (abs_grad_phi+dx*dx*(abs_grad_phi == 0)));
H1_abs = H1_abs(3:end-2,3:end-2);
H2_abs = H2_abs(3:end-2,3:end-2);

delta = Vn.*abs_grad_phi;
delta = delta(3:end-2,3:end-2);










