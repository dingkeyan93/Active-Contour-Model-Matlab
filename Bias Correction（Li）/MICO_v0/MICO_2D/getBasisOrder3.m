% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://imagecomputing.org/~cmli/
function B = getBasisOrder3(Height,Wide)

for i =1:Height
    x(i,:) = -1:2/(Wide-1):1;
end
for i =1:Wide
    temp = -1:2/(Height-1):1;
    y(:,i) = temp';
end


bais = zeros(Height,Wide,10);
bais(:,:,1) = 1;
bais(:,:,2) = x;
bais(:,:,3) = (3.*x.*x - 1)./2;
bais(:,:,4) = (5.*x.*x.*x - 3.*x)./2;
bais(:,:,5) = y;
bais(:,:,6) = x.*y;
bais(:,:,7) = y.*(3.*x.*x -1)./2;
bais(:,:,8) = (3.*y.*y -1)./2;
bais(:,:,9) = (3.*y.*y -1).*x./2;
bais(:,:,10) = (5.*y.*y.*y -3.*y)./2;


B = bais;
for kk=1:10
    A=bais(:,:,kk).^2;
    r = sqrt(sum(A(:)));
    B(:,:,kk)=bais(:,:,kk)/r;
end
