function [bdy,bdx]=backward_gradient(f);
% function [bdx,bdy]=backward_gradient(f);
%  
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li
[nr,nc]=size(f);
bdx=zeros(nr,nc);
bdy=zeros(nr,nc);

bdx(2:nr,:)=f(2:nr,:)-f(1:nr-1,:);
bdy(:,2:nc)=f(:,2:nc)-f(:,1:nc-1);