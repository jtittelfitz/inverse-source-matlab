   
%%
%% fine or coarse mesh?
%%
    h = 3/res; %h = 0.0125;%0.025;
    k = h/5;   %(h/5) k = 0.0025;%0.005;%
    
%%
%% initialize grid
%%
	leftB = -2;             % left and right endpoints
	rightB = 2;             %   of smaller region
    left = -5;              % left and right endpoints
    right = 5;              %   of "Rn"

    t = 0:k:T;                          % time interval
	[x1,x2] = meshgrid(left:h:right);   % "Rn"
    X = invsource_grid(x1,x2,t,h,k,c,scale);
	[y1,y2] = meshgrid(leftB:h:rightB); % smaller region
    Y = invsource_grid(y1,y2,t,h,k,c,scale);
    
    C = speeds(c,x1,x2,scale);
    %c_11 = C{1,1};
    %c_12 = C{1,2};
    %c_13 = C{1,3};
    %c_22 = C{2,2};
    %c_23 = C{2,3};
    %c_33 = C{3,3};
    	
    supp = (x1 >= -1.5) & (x1 < 1.5) & (x2 >= -1.5) & (x2 < 1.5);
    suppv = (y1 >= -1.5) & (y1 < 1.5) & (y2 >= -1.5) & (y2 < 1.5);
    omega = (x1 >= -2) & (x1 <= 2) & (x2 >= -2) & (x2 <= 2);
    omegav = (y1 >= -2) & (y1 <= 2) & (y2 >= -2) & (y2 <= 2);
    
    slice = (y2 == 0);