function [KE,PE,CE] = display_energies( u,u_t,X,t_j,display )
    
    if nargin < 4
        display = 1;
    end

    %if iscell(u)
     %   for j = 1 : length(u)
    
    KE = kinetic_energy(u_t,X);
    PE = zeros;
    %PE = potential_energy(u,X);
    CE = concentration(u,u_t,X);
    if display
        figure
        subplot(2,2,1); plot(X.t,KE); title('Kinetic Energy'); vline(t_j,'r');
        subplot(2,2,2); plot(X.t,CE); title('Volume of Support'); vline(t_j,'r'); %plot(X.t,PE); title('Potential Energy'); vline(t_j,'r');
        subplot(2,2,3); plot(X.t,KE./CE); title('Kinetic Energy Concentration'); vline(t_j,'r'); %plot(X.t,KE + PE); title('Total Energy');
    end
    
    %c_max
    %ke_density = ((0.5*u(:,:,t_max+1)-0.5*u(:,:,t_max-1))/X.k).^2;
    
    %odd_diff = (u(:,:,t_max+2)-0.5*u(:,:,t_max+1))/X.k + (u(:,:,t_max)-0.5*u(:,:,t_max-1))/X.k;
    %max(max(ke_density))
    %subplot(2,2,4); imagesc(   ((u(:,:,t_max)).^2 < (X.k)^2) & ( ((u(:,:,t_max) - u(:,:,t_max-1))/X.k).^2 > X.k)   ) %imagesc(odd_diff)%imagesc(ke_density > c_max - X.k)
    
    
end

