%/* this function computes the signed distance as initial level set function  */
function u = signed_distance(I,xcontour, ycontour,margin)
% I is the image matrix
% nrow is the no of rows
% ncol is the no of columns

[nrow, ncol] = size(I);
[temp, contsize] = size(xcontour);

Mark = zeros(nrow, ncol);

for y=1:nrow,
    for x=1:ncol
        if  (x > ncol-margin+1) | (x < margin) | (y < margin) | (y > nrow-margin+1)
            Mark(y,x) = -1;
        end;
    end;
end;

for y = 1:nrow,
    for x =1: ncol,
        u(y,x) = sqrt(min((x-xcontour).^2+(y-ycontour).^2));
        if Mark(y,x) == -1
            u(y,x) = -u(y,x);
        end;
    end;
end;



