classdef TAT_grid
    %TAT_grid a simple object to keep grid and grid spacing together
    
    properties
        x1
        x2
        t
        h
        k
    end
    
    methods
        function g = TAT_grid(x1,x2,t,h,k)
            g.x1 = x1;
            g.x2 = x2;
            g.t = t;
            g.h = h;
            g.k = k;
        end
    end
    
end

