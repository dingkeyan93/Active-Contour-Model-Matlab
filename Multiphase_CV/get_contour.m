%/* this function computes the phi required initially  */
function [xcontour, ycontour] = get_phi(I, nrow, ncol,margin)
% I is the image matrix
% nrow is the no of rows
% ncol is the no of columns


count=1;
x=margin;
for y=margin:nrow-margin+1,
   xcontour(count) = x;
   ycontour(count) = y;
   count=count+1;
end;
y=nrow-margin+1;
for x=margin+1:ncol-margin+1,
   xcontour(count) = x;
   ycontour(count) = y;
   count=count+1;
end;
     
x=ncol-margin+1;
for y=nrow-margin:-1:margin,
   xcontour(count) = x;
   ycontour(count) = y;
   count=count+1;
end;

y=margin;
for x=ncol-margin:-1:margin+1,
   xcontour(count) = x;
   ycontour(count) = y;
   count=count+1;
end;



