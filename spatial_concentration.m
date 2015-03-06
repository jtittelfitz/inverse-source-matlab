function v = spatial_concentration(u,t,X,r)
    if nargin < 3
        r = 10;
    end
    v = zeros(length(X.x1), length(X.x2));
       
    w = ((u(:,:,floor(t/X.k)) - u(:,:,floor(t/X.k)-1))/X.k).^2;
        
    f = @(x,y) full(sum(sum(w.*sparse(disc(x,y,r,X)))));
    v = arrayfun(f,X.x1,X.x2);
        
    %A = (abs(u(:,:,n)) > tol) | (abs(u(:,:,n) - u(:,:,n-1))/X.k > tol);
    %CE(n) = (X.h)^2*sum(sum(  A  ));    

end

