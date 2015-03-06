function m = mask(u,u_t,t,zero_tol,nonzero_tol)
% mask - finds support of sources
%   finds locations where displacement vanishes, velocity does not

    m = ((u(:,:,t)).^2 < zero_tol) & ( (u_t(:,:,t)).^2 > nonzero_tol);


end

