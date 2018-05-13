% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://imagecomputing.org/~cmli/
function B = getBasisOrder4(Height,Wide)

for i =1:Height
    x(i,:) = -1:2/(Wide-1):1;
end
for i =1:Wide
    temp = -1:2/(Height-1):1;
    y(:,i) = temp';
end

bais = zeros(Height,Wide,15);
bais(:,:,1) = 1;
bais(:,:,2) = x;
bais(:,:,3) = (3.*x.*x - 1)./2;
bais(:,:,4) = (5.*x.*x.*x - 3.*x)./2;
bais(:,:,5) = (35.*x.*x.*x.*x - 30.*x.*x+3)./8;
bais(:,:,6) = y;
bais(:,:,7) = x.*y;
bais(:,:,8) = (3.*x.*x - 1).*y./2;
bais(:,:,9) = (5.*x.*x.*x - 3.*x).*y./2;
bais(:,:,10) = (3.*y.*y - 1)./2;
bais(:,:,11) = (3.*y.*y - 1).*x./2;
bais(:,:,12) = (3.*x.*x - 1).*(3.*y.*y - 1)./4;
bais(:,:,13) = (5.*y.*y.*y- 3.*y)./2;
bais(:,:,14) = (5.*y.*y.*y- 3.*y).*x./2;
bais(:,:,15) = (35.*y.*y.*y.*y - 30.*y.*y+3)./8;




for kk=1:15
    A=bais(:,:,kk).^2;
    r = sqrt(sum(A(:)));
    B(:,:,kk)=bais(:,:,kk)/r;
end
