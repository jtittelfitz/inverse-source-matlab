function u_t = velocity(u,X,diff_type)
    if nargin < 3
        diff_type = 'center';
    end

    if strcmp(diff_type,'forward')
        u_t = (u - circshift(u,[0,0,1]))/X.k;
    elseif strcmp(diff_type,'backward')
        u_t = (circshift(u,[0,0,-1]-u))/X.k;
    else % use centered difference
        u_t = (circshift(u,[0,0,-1]) - circshift(u,[0,0,1]))/(2*X.k);
    end

    u_t(:,:,1) = zeros;
    u_t(:,:,length(X.t)) = zeros;
end

