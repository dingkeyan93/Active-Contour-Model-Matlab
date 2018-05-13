function [M_out, C_out]=sortMemC(M, C)

[C_out IDX]=sort(C);

for k = 1 : length(C)
    M_out(:,:,:,k) = M(:,:,:,IDX(k));
end

