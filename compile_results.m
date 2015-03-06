results = zeros(length(separation_and_energy),3);
for test = 1 : length(separation_and_energy)
    separation = separation_and_energy(test,1);
    results(test,2) = separation;
    rel_energy = separation_and_energy(test,2);    
    results(test,3) = rel_energy;
    if ((rel_energy < 0.01) & separation > 0) | ((rel_energy >= 0.01) & separation <= 0)
        results(test,1) = 1;
    end
end
       