function u = re_init_u( u )
% re-initialize the u funcion

% get the  size of u
[M,N] = size (u);

[x,y] = find( isfront( u ) );
L = length ( x );
for i = 1 : M
    for j = 1 : N
        min_dist = sqrt(M^2+N^2);
        for k= 1 : L
            dist = sqrt((i-x(k))^2+(j-y(k))^2) ;
            if dist<min_dist
                min_dist = dist;
            end
        end
        if u(i,j) > 0 
            u(i,j) = min_dist;
        else
            u(i,j) = - min_dist;
        end
    end
end


