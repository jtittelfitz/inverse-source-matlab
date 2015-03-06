function KE = kinetic_energy(u_t,X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    KE = zeros(1,length(X.t));
    KE(1) = 0;    
    for n = 2 : length(X.t) -1
        KE(n) = (X.h)^2*norm(u_t(:,:,n),'fro')^2;
    end
    KE(length(X.t)) = 0;
end

