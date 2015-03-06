%%
%% Key Program Parameters
%%

res = 120;       % resolution of initial data
c = 2;          % choice of wave speed (see speeds.m)
scale = 1;      % scale
T = 2*scale;    % maximum value of T
video = true;       % show movies of wave propagation

%%
%% initialize grid
%%   

gridspeedinit;    
diffs = diff_scheme_obj(h,length(x1));


%% initialize solution
u = zeros(length(x1),length(x2),length(t));

% f and g are initial displacement and velocity for corresponding initial value problems (via Duhamel's)
f = zeros(length(x1),length(x2));
g = double(disc([0,0],0.5,X));            

% solve wave equation forward in time
u = forward(f,g,X,diffs);

if video; play(u,velocity(u,X),k,1,0.1); end        

toc
