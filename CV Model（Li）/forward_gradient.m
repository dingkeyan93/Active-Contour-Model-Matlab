function [fdy,fdx]=forward_gradient(f);
% function [fdx,fdy]=forward_gradient(f);
%  
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li
[nr,nc]=size(f);
fdx=zeros(nr,nc);
fdy=zeros(nr,nc);

a=f(2:nr,:)-f(1:nr-1,:);
fdx(1:nr-1,:)=a;
b=f(:,2:nc)-f(:,1:nc-1);
fdy(:,1:nc-1)=b;
