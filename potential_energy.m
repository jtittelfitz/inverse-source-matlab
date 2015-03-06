function PE = potential_energy(u,X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    N_grid = length(X.x1);
    [DX,DZ] = quick_diff_scheme( X.h,N_grid);
    PE = zeros(1,length(X.t));
    for n = 1 : length(X.t)
        v = reshape(u(:,:,n),N_grid^2,1);
        PE(n) = (X.h)^2*(((X.c).^2).*(norm(DX*v,'fro')^2 + norm(DZ*v,'fro')^2));
    end
end

