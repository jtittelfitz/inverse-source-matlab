function v = reverse (phi,Lambda_f,C,X,Y,record)
    
    if nargin < 6
        record = 1;
    end

    phi_1 = phi{1};
    phi_2 = phi{2};
    
    Lambda_f_u1 = Lambda_f{1};
    Lambda_f_u2 = Lambda_f{2};
    
    y1 = Y.x1;
    y2 = Y.x2;
    t = Y.t;
    h = Y.h;
    k = Y.k;
    
    x1 = X.x1;
    x2 = X.x2;
    omega = (x1 >= -2) & (x1 <= 2) & (x2 >= -2) & (x2 <= 2);
    
    N_grid = length(y1);
    
    c_11_0 = reshape(C{1,1}(omega),N_grid^2,1);
    c_22_0 = reshape(C{2,2}(omega),N_grid^2,1);
    c_12_0 = reshape(C{1,2}(omega),N_grid^2,1);
    c_33_0 = reshape(C{3,3}(omega),N_grid^2,1);
    
    [DX,DZ,DRX,DRZ,DLX,DLZ,DPHX,DPHZ,DMHX,DMHZ] = diff_scheme(h,N_grid);

    % initialize solution array v

    v_1_1 = zeros(N_grid^2,1);
    v_1_2 = zeros(N_grid^2,1);
    v_1_3 = zeros(N_grid^2,1);
    v_2_1 = zeros(N_grid^2,1);
    v_2_2 = zeros(N_grid^2,1);
    v_2_3 = zeros(N_grid^2,1);

    % v(:,:,1) is previous time step (left = 0)
    % v(:,:,2) is current time step
    v_1_2 = reshape(phi_1,N_grid^2,1);                   
    v_2_2 = reshape(phi_2,N_grid^2,1);

    bdyv = (((abs(y1) == 2) & (abs(y2)<=2)) | ((abs(y2) == 2) & (abs(y1) <=2)));

    if record == 2
        record_v_1 = zeros(length(y1),length(y2),length(t));                                    
        record_v_1(:,:,1) = phi_1;
        record_v_2 = zeros(length(y1),length(y2),length(t));                                    
        record_v_2(:,:,1) = phi_2;
    end

    clear phi_1;
    clear phi_2;

    % v(:,:,3) is next time step

    v_1_3 = v_1_2 + 0.5*k^2*((DPHX*(c_11_0)).*(DRX*v_1_2)...
        - (DMHX*(c_11_0)).*(DLX*v_1_2)...
        + (DPHZ*c_33_0).*(DRZ*v_1_2) - (DMHZ*c_33_0).*(DLZ*v_1_2)...
        + DX*(c_12_0.*(DZ*v_2_2)) + DZ*(c_33_0.*(DX*v_2_2)));
    v_2_3 = v_2_2 + 0.5*k^2*((DPHZ*(c_22_0)).*(DRZ*v_2_2)...
        - (DMHZ*(c_22_0)).*(DLZ*v_2_2)...
        + (DPHX*c_33_0).*(DRX*v_2_2) - (DMHX*c_33_0).*(DLX*v_2_2)...
        + DZ*(c_12_0.*(DX*v_1_2)) + DX*(c_33_0.*(DZ*v_1_2)));


    v_1 = reshape(v_1_3,N_grid,N_grid);
    v_2 = reshape(v_2_3,N_grid,N_grid);
    v_1(bdyv) = Lambda_f_u1(:,length(t) - 1);
    v_2(bdyv) = Lambda_f_u2(:,length(t) - 1);

    v_1_3 = reshape(v_1,N_grid^2,1);
    v_2_3 = reshape(v_2,N_grid^2,1);

    % update arrays in preparation for first pass
    v_1_1 = v_1_2;
    v_1_2 = v_1_3;
    v_2_1 = v_2_2;
    v_2_2 = v_2_3;

    if record == 2
        record_v_1(:,:,2) = reshape(v_1_2,N_grid,N_grid);
        record_v_2(:,:,2) = reshape(v_2_2,N_grid,N_grid);
    end

    for n = 2: length(t)-1
        v_1_3 = 2*v_1_2 - v_1_1 + k^2*((DPHX*(c_11_0)).*(DRX*v_1_2)...
            - (DMHX*(c_11_0)).*(DLX*v_1_2)...
            + (DPHZ*c_33_0).*(DRZ*v_1_2) - (DMHZ*c_33_0).*(DLZ*v_1_2)...
            + DX*(c_12_0.*(DZ*v_2_2)) + DZ*(c_33_0.*(DX*v_2_2)));
        v_2_3 = 2*v_2_2 - v_2_1 + k^2*((DPHZ*(c_22_0)).*(DRZ*v_2_2)...
            - (DMHZ*(c_22_0)).*(DLZ*v_2_2)...
            + (DPHX*c_33_0).*(DRX*v_2_2) - (DMHX*c_33_0).*(DLX*v_2_2)...
            + DZ*(c_12_0.*(DX*v_1_2)) + DX*(c_33_0.*(DZ*v_1_2)));

        v_1 = reshape(v_1_3,N_grid,N_grid);
        v_2 = reshape(v_2_3,N_grid,N_grid);
        v_1(bdyv) = Lambda_f_u1(:,length(t) - n);
        v_2(bdyv) = Lambda_f_u2(:,length(t) - n);

        v_1_3 = reshape(v_1,N_grid^2,1);
        v_2_3 = reshape(v_2,N_grid^2,1);

        if record == 1
            imagesc(v_1); drawnow
        end
        % update arrays for next pass
        v_1_1 = v_1_2;
        v_1_2 = v_1_3;
        v_2_1 = v_2_2;
        v_2_2 = v_2_3;
    end
    v_1 = reshape(v_1_3,N_grid,N_grid);
    v_2 = reshape(v_2_3,N_grid,N_grid);

    v = {v_1,v_2};
end
