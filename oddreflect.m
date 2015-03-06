function w = oddreflect(u,t_0,T,X)
    w = zeros(length(X.x1),length(X.x2),floor(T/X.k)+1);
    start = floor(t_0/X.k) + 1;
    finish = floor(T/X.k) + 1;
    w(:,:,start:finish) = u(:,:,1:finish-start+1);
    temp = flipdim(u,3);
    w(:,:,1:start - 1) = -temp(:,:,length(X.t) - start + 1:length(X.t)-1);    
end