classdef invsource_grid
    %TAT_grid a simple object to keep grid and grid spacing together
    
    properties
        x1
        x2
        t
        h
        k
        C
    end
    
    methods
        function g = invsource_grid(x1,x2,t,h,k,c,scale)
            g.x1 = x1;
            g.x2 = x2;
            g.t = t;
            g.h = h;
            g.k = k;
            g.C = speeds(c,x1,x2,scale);
        end
    end
    
end

