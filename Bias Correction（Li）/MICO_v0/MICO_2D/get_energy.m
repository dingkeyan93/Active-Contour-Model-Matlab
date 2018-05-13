% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://imagecomputing.org/~cmli/
function energy = get_energy(Img,b,C,M,ROI,q)

N = size(M,3);
energy = 0;

for k = 1 : N
    C_k=C(k)*ones(size(Img));
    energy = energy + sum(sum((Img.*ROI - b.*C_k.*ROI).^2.*M(:,:,k).^q));
end