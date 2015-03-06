function hhh = component_energies(w,X,energy_type)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    figure
    colors = ['r','b','g'];
    for j = 1 : length(w)
        w_t = velocity(w{j},X);
        KE = kinetic_energy(w_t,X);
        CE = concentration(w{j},w_t,X);
        k = mod(j,length(colors)) + 1;
        if strcmp(energy_type,'kinetic')
            plot(X.t,KE,colors(k))    
            hold
        else
            plot(X.t,KE./CE,colors(k))
            hold
        end
    end
    hold off
end

