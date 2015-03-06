function CE = concentration(u,u_t,X,tol)
    if nargin < 4
        tol = X.k;
    end
    CE = zeros(1,length(X.t));
    CE(1) = 1;
    for n = 2 : length(X.t) -1
        A = (abs(u(:,:,n)) > tol) | (abs(u_t(:,:,n)) > tol);
        CE(n) = (X.h)^2*sum(sum(  A  ));
    end
    CE(length(X.t)) = 1;
end

