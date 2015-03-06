
tic
%%
%% Key Program Parameters
%%
	
	res = 120;       % resolution of initial data
    c = 1;          % wave speed
    scale = 1;      % scale
	T = 2*scale;    % maximum value of T
    N = 1;          % number of sources
    
    zero_tol = 1.0e-10; % criteria for determining support of a source
    nonzero_tol = 0.5;
    mask_tol = 1;
    mask_search_width = 2;
    
    video = true;       % show movies of wave propagation
    display = true;    % show graphs of energies
    output = true;     % text output
    write = false;      % catalog results
        test_energy = true; % catalog which type of failures?
        test_mask = true;
    randomize_input = false;
    
    num_tests = 1;   
    results = zeros(num_tests,1);    separation_and_energy = zeros(num_tests,2);

%%
%% initialize output file
%%
if write
    path = 'output/';
    mkdir(path);
    file = sprintf('%s/output.txt',path);
    output_file = fopen(file,'w');			
end
    
%%
%% initialize grid
%%    
gridspeedinit;    
diffs = diff_scheme_obj(h,length(x1));
       

    
    %%	
    %% set (randomize) initial data
    %%    
    mask_flag = false;
    u = zeros(length(x1),length(x2),length(t)); 
    
    f = zeros(length(x1),length(x2));%double(disc(centers(j,:),r(j),X));            %zeros(length(x1),length(x2));                 
    g = sin(6*x1).*sin(4*x2).*(abs(x1) < pi/2 & abs(x2) < pi/2);%double(disc([0,0],1,X));            

    u = forward(f,g,X,diffs);
    u_t = velocity(u,X);
    play(u,k,0.1,u_t);        
    



toc