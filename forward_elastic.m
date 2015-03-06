function Lambda_f = forward (f,C,X,mesh,Z,record,prec,errortest)     

    if nargin < 9
        record = 1;
        prec = 1;
        errortest = 0;
    end
    
    x1 = X.x1;
    x2 = X.x2;
    t = X.t;
    h = X.h;
    k = X.k;
    
% initialize spatial derivative operators
    N_grid = length(x1);
    
	if prec			
		[DX_11_DX,DZ_33_DZ,DX_12_DZ,DZ_33_DX,DZ_22_DZ,DX_33_DX,DZ_12_DX,DX_33_DZ] = diff_scheme_prec( C, h, N_grid);
	else
		[DX,DZ,DRX,DRZ,DLX,DLZ,DPHX,DPHZ,DMHX,DMHZ] = diff_scheme( h,N_grid );
		c_11 = reshape(C{1,1},N_grid^2,1);
		c_12 = reshape(C{1,2},N_grid^2,1);
		c_13 = reshape(C{1,3},N_grid^2,1);
		c_22 = reshape(C{2,2},N_grid^2,1);
		c_23 = reshape(C{2,3},N_grid^2,1);
		c_33 = reshape(C{3,3},N_grid^2,1);
	end

    f_1 = reshape(f{1},N_grid^2,1);
    f_2 = reshape(f{2},N_grid^2,1);

    % setup to calculate solution using 2nd order method

    u_1_1 = zeros(N_grid^2,1);
    u_1_2 = f_1;

    u_2_1 = zeros(N_grid^2,1);
    u_2_2 = f_2;
    if prec ==1
    u_1_3 = u_1_2 + 0.5*k^2*(DX_11_DX*u_1_2 + DZ_33_DZ*u_1_2...
            + DX_12_DZ*u_2_2 + DZ_33_DX*u_2_2);

    u_2_3 = u_2_2 + 0.5*k^2*(DZ_22_DZ*u_2_2 + DX_33_DX*u_2_2...
            + DZ_12_DX*u_1_2 + DX_33_DZ*u_1_2);

    else
    u_1_3 = u_1_2 + 0.5*k^2*((DPHX*(c_11)).*(DRX*u_1_2)...
        - (DMHX*(c_11)).*(DLX*u_1_2)...
        + (DPHZ*c_33).*(DRZ*u_1_2) - (DMHZ*c_33).*(DLZ*u_1_2)...
        + DX*(c_12.*(DZ*u_2_2)) + DZ*(c_33.*(DX*u_2_2)));

    u_2_3 = u_2_2 + 0.5*k^2*((DPHZ*(c_22)).*(DRZ*u_2_2)...
        - (DMHZ*(c_22)).*(DLZ*u_2_2)...
        + (DPHX*c_33).*(DRX*u_2_2) - (DMHX*c_33).*(DLX*u_2_2)...
        + DZ*(c_12.*(DX*u_1_2)) + DX*(c_33.*(DZ*u_1_2)));
    end

    u_1_1 = u_1_2;
    u_1_2 = u_1_3;
    u_2_1 = u_2_2;
    u_2_2 = u_2_3;

    f_1 = reshape(f_1,N_grid,N_grid);
    f_2 = reshape(f_2,N_grid,N_grid);

    %initialize boundary data

    if mesh == 1
        z1 = Z.x1;
        z2 = Z.x2;
        bdy = (((abs(z1) == 2) & (abs(z2)<=2)) | ((abs(z2) == 2) & (abs(z1) <=2)));
        bdyx = (((abs(x1) == 2) & (abs(x2)<=2)) | ((abs(x2) == 2) & (abs(x1) <=2)));
        Lambda_f_u1 = zeros(length(z1(bdy)),length(t));
        Lambda_f_u2 = zeros(length(z1(bdy)),length(t));
    else
        bdy = (((abs(x1) == 2) & (abs(x2)<=2)) | ((abs(x2) == 2) & (abs(x1) <=2)));
        Lambda_f_u1 = zeros(length(x1(bdy)),length(t));
        Lambda_f_u2 = zeros(length(x1(bdy)),length(t));
    end

    % Lambda_f(x,y,n) = u(x,y,(n-1)k) for x,y in boundary of Omega                                     
    u_1 = reshape(u_1_3,N_grid,N_grid);
    u_2 = reshape(u_2_3,N_grid,N_grid);

    if mesh == 1
        g_1 = downsample(downsample(f_1,2)',2)';
        g_2 = downsample(downsample(f_2,2)',2)';
        h_1 = downsample(downsample(u_1,2)',2)';
        h_2 = downsample(downsample(u_2,2)',2)';    
        Lambda_f_u1(:,1) = g_1(bdy);
        Lambda_f_u2(:,1) = g_2(bdy);
        Lambda_f_u1(:,2) = h_1(bdy);
        Lambda_f_u2(:,2) = h_2(bdy);
        clear g_1;
        clear g_2;
        clear h_1;
        clear h_2;
    else
        Lambda_f_u1(:,1) = f_1(bdy);
        Lambda_f_u2(:,1) = f_2(bdy);
        Lambda_f_u1(:,2) = u_1(bdy);
        Lambda_f_u2(:,2) = u_2(bdy);   
    end

    if record == 1
        %mov = avifile('mesh2.avi');
        %colormap(gray)
        %if mesh == 1
        %    %imagesc(u_1(omega))%(152:161,152:161,n))%
        %    subplot(1,2,1)
        %    plot3(x1(bdyx),x2(bdyx),u_1(bdyx),'.')
        %    subplot(1,2,2)
        %    plot3(z1(bdy),z2(bdy),Lambda_f_u1(:,1),'.')
        %    %caxis([-1,1])
        %    %colorbar
        %else
        %    plot3(x1(bdy),x2(bdy),u_1(bdy),'.')
        %end    
	imagesc(u_1)
        drawnow
        %frame = getframe(gcf);
        %mov = addframe(mov,frame);
        %record_1 = zeros(length(x1),length(x2),length(t));
        %record_2 = zeros(length(x1),length(x2),length(t));
        %record_1(:,:,1) = f_1;
        %record_2(:,:,1) = f_2;
        %record_1(:,:,2) = u_1;
        %record_2(:,:,2) = u_2;
    end


    for n = 2: length(t)-1
        %finite difference scheme

        %% 2nd Order Method
        if prec ==1 
        u_1_3 = 2*u_1_2 - u_1_1 + k^2*(DX_11_DX*u_1_2 + DZ_33_DZ*u_1_2...
            + DX_12_DZ*u_2_2 + DZ_33_DX*u_2_2);

        u_2_3 = 2*u_2_2 - u_2_1 + k^2*(DZ_22_DZ*u_2_2 + DX_33_DX*u_2_2...
            + DZ_12_DX*u_1_2 + DX_33_DZ*u_1_2);
        else
        u_1_3 = 2*u_1_2 - u_1_1 + k^2*((DPHX*(c_11)).*(DRX*u_1_2)...
            - (DMHX*(c_11)).*(DLX*u_1_2)...
            + (DPHZ*c_33).*(DRZ*u_1_2) - (DMHZ*c_33).*(DLZ*u_1_2)...
            + DX*(c_12.*(DZ*u_2_2)) + DZ*(c_33.*(DX*u_2_2)));
        u_2_3 = 2*u_2_2 - u_2_1 + k^2*((DPHZ*(c_22)).*(DRZ*u_2_2)...
            - (DMHZ*(c_22)).*(DLZ*u_2_2)...
            + (DPHX*c_33).*(DRX*u_2_2) - (DMHX*c_33).*(DLX*u_2_2)...
            + DZ*(c_12.*(DX*u_1_2)) + DX*(c_33.*(DZ*u_1_2)));
        end
        %store boundary data
        u_1 = reshape(u_1_3,N_grid,N_grid);
        u_2 = reshape(u_2_3,N_grid,N_grid);

        if mesh == 1
            g_1 = downsample(downsample(u_1,2)',2)';
            g_2 = downsample(downsample(u_2,2)',2)';
            Lambda_f_u1(:,n+1) = g_1(bdy);
            Lambda_f_u2(:,n+1) = g_2(bdy);
            clear g_1;
            clear g_2;
        else
            Lambda_f_u1(:,n+1) = u_1(bdy);
            Lambda_f_u2(:,n+1) = u_2(bdy);
        end

        if record == 1
            %imagesc(u_1(omega))%(152:161,152:161,n))%
            %if mesh == 1
            %    subplot(1,2,1)
            %    plot3(x1(bdyx),x2(bdyx),u_1(bdyx),'.')
            %    subplot(1,2,2)
            %    plot3(z1(bdy),z2(bdy),Lambda_f_u1(:,n+1),'.')
            %    caxis([-1,1])
            %else
            %    plot3(x1(bdy),x2(bdy),u_1(bdy),'.')
            %end    
            %colorbar
		imagesc(u_1)
            drawnow
            %frame = getframe(gcf);
            %mov = addframe(mov,frame);
            %record_1(:,:,n+1) = u_1;
            %record_2(:,:,n+1) = u_2;	
        end

        % update arrays in preparation for first pass
        u_1_1 = u_1_2;
        u_1_2 = u_1_3;
        u_2_1 = u_2_2;
        u_2_2 = u_2_3;
    end

    if mesh == 1
        for n = 2 : floor(length(t)/2) + 1
            Lambda_f_u1(:,n) =  Lambda_f_u1(:,2*n - 1);
            Lambda_f_u2(:,n) =  Lambda_f_u2(:,2*n - 1);
        end
    end
    
    Lambda_f = {Lambda_f_u1,Lambda_f_u2};
    
    if record == 1
        %mov = close(mov);
    end
    
    if errortest
        if mesh 
            u_1 = downsample(downsample(u_1,2)',2)';
            u_2 = downsample(downsample(u_2,2)',2)';
        end
    else
        clear u_1;
        clear u_2;
    end

end
