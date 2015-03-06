function [ C ] = speeds( c,x1,x2,scale,elastic )
%speeds returns preset lame parameters
% c_12 = "lambda"
% c_33 = "mu"

    if (nargin < 5)
        scale = 1;
        elastic = false;
    end
    
    if ~elastic
        C = zeros(length(x1),length(x2));

        %% homogeneous;
        if c == 1
            C = (1/scale^2)*1*ones(length(x1),length(x2));            
        end

        %% non-homogeneous;
        if c == 2
        %% current: 
            C = (1/scale^2)*(1.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2; %(1/scale^2)*(1.0 + 0.02*sin(2*pi*x1) + 0.01*cos(2*pi*x2)).^2;
            %C = (1/scale^2)*(2.0 + 0.1*cos(2*pi*x1) + 0.2*sin(2*pi*x2)).^2;

            cutoff_value = (1/scale^2)*1;            
            C = cutoff(C,x1,x2,cutoff_value);
            
        end	

        %% c_P > c_S; (1,4); T_0 = 1.8068
        if c == 3            
            C = (1/scale^2)*(1.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2;            
            cutoff_value = (1/scale^2)*1;            
            C = cutoff(C,x1,x2,cutoff_value);
        end

        %% trapping speed; (2,3); T_0 = 1.7
        if c == 4            
            C = (1/scale^2)*(1.25 + sin(2.0*pi*x1).*cos(2*pi*x2)).^2;            
            cutoff_value = (1/scale^2)*1.25^2;            
            C = cutoff(C{3,3},x1,x2,cutoff_value);
        end	
    else
        C = cell(3,3);    

        %% homogeneous;
        if c == 1
            C{1,2} = (1/scale^2)*1*ones(length(x1),length(x2));
            C{3,3} = zeros(length(x1),length(x2));%(1/scale^2)*4*ones(length(x1),length(x2));
        end

        %% non-homogeneous;
        if c == 2
        %% current: 
            C{1,2} = (1/scale^2)*(1.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2;
            C{3,3} = (1/scale^2)*(2.0 + 0.1*cos(2*pi*x1) + 0.2*sin(2*pi*x2)).^2;

            %% experimental:
            %c_12 = (1/scale^2)*(2.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2;
            %c_33 = (1/scale^2)*(2.0 + 0.1*cos(2*pi*x1) + 0.2*sin(2*pi*x2)).^2;

            %c_12 = (1/scale^2)*(2.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2;
            %c_33 = (1.0 + 0.1*cos(2*pi*x1) + 0.2*sin(2*pi*x2)).^2;

            cutoff_value_12 = (1/scale^2)*1;
            cutoff_value_33 = (1/scale^2)*4;
            C{1,2} = cutoff(C{1,2},x1,x2,cutoff_value_12);
            C{3,3} = cutoff(C{3,3},x1,x2,cutoff_value_33);
        end	

        %% c_P > c_S; (1,4); T_0 = 1.8068
        if c == 3
            C{1,2} = (1/scale^2)*14*ones(length(x1),length(x2)); 
            C{3,3} = (1/scale^2)*(1.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2;
            cutoff_value_12 = (1/scale^2)*14;
            cutoff_value_33 = (1/scale^2)*1;
            C{1,2} = cutoff(C{1,2},x1,x2,cutoff_value_12);
            C{3,3} = cutoff(C{3,3},x1,x2,cutoff_value_33);
        end

        %% trapping speed; (2,3); T_0 = 1.7
        if c == 4
            C{1,2} = (1/scale^2)*1*ones(length(x1),length(x2));
            C{3,3} = (1/scale^2)*(1.25 + sin(2.0*pi*x1).*cos(2*pi*x2)).^2;
            cutoff_value_12 = (1/scale^2)*1;
            cutoff_value_33 = (1/scale^2)*1.25^2;
            C{1,2} = cutoff(C{1,2},x1,x2,cutoff_value_12);
            C{3,3} = cutoff(C{3,3},x1,x2,cutoff_value_33);
        end	

        %% c^+ >> c^-
        if c == 5
            C{1,2} = (1/scale^2)*(2.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2;
            C{3,3} = (1/scale^2)*(0.001)^2*ones(length(x1),length(x2)); 
            cutoff_value_12 = (1/scale^2)*4;
            cutoff_value_33 = (1/scale^2)*(0.001)^2;
            C{1,2} = cutoff(C{1,2},x1,x2,cutoff_value_12);
            C{3,3} = cutoff(C{3,3},x1,x2,cutoff_value_33);
        end

        % isotropic, c_11 = c_22 = c_12 + 2*c_33
        if c < 6
            C{1,1} = C{1,2} + 2*C{3,3};
            C{2,2} = C{1,2} + 2*C{3,3};
        elseif c == 6
            C{1,2} = (1/scale^2)*(1.0 + 0.2*sin(2*pi*x1) + 0.1*cos(2*pi*x2)).^2;
            C{3,3} = (1/scale^2)*(2.0 + 0.1*cos(2*pi*x1) + 0.2*sin(2*pi*x2)).^2;
            cutoff_value_12 = (1/scale^2)*1;
            cutoff_value_33 = (1/scale^2)*4;
            C{1,2} = cutoff(C{1,2},x1,x2,cutoff_value_12);
            C{3,3} = cutoff(C{3,3},x1,x2,cutoff_value_33);
            C{1,1} = C{1,2} + 2*C{3,3} + 0.25;
            C{2,2} = C{1,2} + 2*C{3,3} - 0.25;
        end

        C{1,3} = sparse(zeros(length(x1),length(x2)));
        C{2,3} = sparse(zeros(length(x1),length(x2)));
    end
end

