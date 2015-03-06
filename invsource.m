
tic
%%
%% Key Program Parameters
%%
	
    res = 120;       % resolution of initial data
    c = 2;          % choice of wave speed (see speeds.m)
    scale = 1;      % scale
     T = 2*scale;    % maximum value of T
    N = 3;          % number of sources
    
    zero_tol = 1.0e-10; % criteria for determining support of a source
    nonzero_tol = 0.5;
    mask_tol = 1;
    mask_search_width = 2;
    
    video = true;       % show movies of wave propagation
    display = false;    % show graphs of energies
    output = true;     % text output
    write = false;      % catalog results to file
        test_energy = true; % catalog which type of failures?
        test_mask = true;
    randomize_input = false; % randomized sources
    
    num_tests = 1;       % parameters for batteries of tests
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
       
for test = 1:num_tests 
    
    %%	
    %% set (randomize) initial data
    %%    t_j's are times of sources
    %%    spatially, sources are characteristic functions of circles, 
    %%    parameterized by center and radius
    
    
     
    
    if randomize_input
        t_j = randi([5,length(X.t)-1-5],1,N)*X.k;
        r = randi(floor(1/X.h),[1,N])*X.h;        
        centers = randi([-ceil(length(Y.x1)/2),floor(length(Y.x2)/2)],N,2)*X.h;        %[-0.28,-1.6];
    else    
        %% presets
        t_j = [1.6850, 1.1150, 1.7850];
        centers = [-1, -1; 1.7250,   -0.3000;   -0.1500,    0.8000];
        r = [0.4250, 0.3750, 0.5000];                
        N = length(t_j);
    end
    
    % initialize solution variables
    % the w's are the waves corresponding to a single source
    % u is the full waveform (all sources)    
    
    w = cell(N);
    solution = zeros(length(x1),length(x2));
    mask_flag = false;
    u = zeros(length(x1),length(x2),length(t));

    % determine space-time separation of sources; if separation is sufficient, proceed
    separation = min(min(stdistances(t_j,centers,r,output)));    
    if separation > 0

        %%
        %% simulate data
        %%
        for j = 1:N
            % f and g are initial displacement and velocity for corresponding initial value problems (via Duhamel's)
            f = zeros(length(x1),length(x2));
            g = double(disc(centers(j,:),r(j),X));            

            % solve wave equation forward in time, then "fake" time-reversal by using Duhamel's/reflection   
            w{j} = forward(f,g,X,diffs);
            w{j} = oddreflect(w{j},t_j(j),T,X);
            
            % add contribution of this single source to full wave
            u = u+w{j};
        end

        % uncomment to save memory; w is not needed from this point  
        %clear w;          
          
        %% visualization/debugging tools

        
        if video; play(u,velocity(u,X),k,1,0.1); end        
        %component_energies(w,X,'kinetic');    
        u_copy = u(:,:,401);   


        %%
        %% reconstruction
        %%    
        for j = 1:N

            %% compute kinetic energy and concentration, find maximum; 
            %%  there might be a source occuring at this time
            u_t = velocity(u,X);

            temp = 1:length(X.t);
            temp_mask = reshape(sum(sum(support_mask(u,u_t,temp,zero_tol,nonzero_tol))),1,length(temp));
            mask_zero = sum(temp_mask == 0);
            fprintf('mask is ident zero at %d of %d times \r',mask_zero,length(temp));
            if (j == 1) && ((length(temp) - mask_zero) > (N - j + 1)); mask_flag = true; disp('bad!'); end            
            %figure; plot(X.t,temp_mask);

            [KE,~,CE] = display_energies(u,u_t,X,t_j,display);  
            if j == 1; orig_energy = max(max(KE)); end

            [c_max,t_max] = max(KE);%max(KE./CE);
            if output; fprintf('Energy concen. max found at %d \r',t_max); end 

            %% find support of the source at this time
            mask = support_mask(u,u_t,t_max,zero_tol,nonzero_tol); %((u(:,:,t_max)).^2 < zero_tol) & ( (u_t(:,:,t_max)).^2 > nonzero_tol);        
            if output; fprintf('Support size: %d \r',sum(sum(mask))); end

            %% check to see if adjacent time works better
            if (sum(sum(mask)) < mask_tol)
                temp = t_max - mask_search_width:t_max + mask_search_width;
                temp(temp < 1) = 1;
                temp(temp > length(X.t)) = length(X.t);
                masks = reshape(sum(sum(support_mask(u,u_t,temp,zero_tol,nonzero_tol))),1,2*mask_search_width+1);

                [a,b] = max(masks);
                t_max = temp(b);
                mask = support_mask(u,u_t,t_max,zero_tol,nonzero_tol);
                if output; fprintf('Support too small; adjusting: t_max =%d , support size =%d \r',t_max,sum(sum(mask))); end
            end        
            if display; subplot(2,2,4); imagesc(mask); end

            %% solve initial value problem related to source estimated above
            g = u_t(:,:,t_max).*mask;
            solution = solution + ((t_max - 1)*X.k)*g;  % not truly a solution; sources in "solution" will get "color-coded" by time they occur
            v = forward(f,g,X,diffs);
            v = oddreflect(v,(t_max-1)*k,T,X);

            %% subtract the effect of this source
            u = u - v;
            play(u,velocity(u,X),k,j+1,0.1);

        end

        %% display computed solution
        if display; figure; imagesc(solution); colorbar; end
        
        %% compute residual (relative) energy to see how well reconstruction worked
        u_t = velocity(u,X);
        [KE,~,CE] = display_energies(u,u_t,X,t_j,display);  
        rel_energy = max(max(KE))/orig_energy;
        separation_and_energy(test,:) = [separation,rel_energy];
        fprintf('Residual (relative) energy: %2.4f \r',rel_energy);

        %% if energy hypotheses were met, but reconstruction failed (or vice versa), document it    
        energy_flag = ((rel_energy >= 0.05) && (separation > 0));
        if write && ((test_energy && energy_flag) || (test_mask && mask_flag)) 
            disp('recorded!');
            fprintf(output_file,'\r -------------------- \r');
            fprintf(output_file,'separation: %2.4f \r',separation);			
            fprintf(output_file,'residual energy (relative): = %2.6f \r',rel_energy);
            fprintf(output_file,'times:\r');
            fprintf(output_file,'%f,\t',t_j);
            fprintf(output_file,'\r\r centers:\r');
            fprintf(output_file,'%f,\t%f;\n',centers');
            fprintf(output_file,'\r radii:\r');
            fprintf(output_file,'%f,\t',r);        
        end
    end
end

if write
    fclose(output_file);
end
toc
