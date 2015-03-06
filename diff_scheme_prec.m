function [DX_11_DX] = diff_scheme_prec( C,h,N_grid ) %,DZ_33_DZ,DX_12_DZ,DZ_33_DX,DZ_22_DZ,DX_33_DX,DZ_12_DX,DX_33_DZ
%diff_scheme instantiate finite difference scheme

	I = speye(N_grid);

    e = ones(N_grid - 1,1);
    D = sparse(diag(e,1));
    D = D - D';
    D(1,1) = -3;
    D(1,2) = 4;
    D(1,3) = -1;
    D(N_grid, N_grid - 2) = 1;
    D(N_grid, N_grid - 1) = -4;
    D(N_grid, N_grid) = 3;
	D = (1/(2*h))*D;

    e = ones(N_grid - 1,1);
    DR = (sparse(diag(e,1)) - I)/(h^2);
    DL = (I - sparse(diag(e,-1)))/(h^2);

    e = ones(N_grid - 1,1);
    DPH = 0.5*(sparse(diag(e,1)) + I);
    DMH = 0.5*(sparse(diag(e,-1)) + I);

    DPHX = kron(DPH,I);
    DPHZ = kron(I,DPH);
    DMHX = kron(DMH,I);
    DMHZ = kron(I,DMH);
    clear DPH;
    clear DMH;

    c_11 = C; %reshape(C{1,1},N_grid^2,1);
    %c_12 = reshape(C{1,2},N_grid^2,1);
    %c_13 = reshape(C{1,3},N_grid^2,1);
    %c_22 = reshape(C{2,2},N_grid^2,1);
    %c_23 = reshape(C{2,3},N_grid^2,1);
    %c_33 = reshape(C{3,3},N_grid^2,1);


        DX_11_DX = spdiags(DPHX*c_11(:),0,N_grid^2,N_grid^2)*kron(DR,I) - spdiags(DPHX*c_11(:),0,N_grid^2,N_grid^2)*kron(DL,I);
        %DZ_33_DZ = spdiags(DPHZ*c_33(:),0,N_grid^2,N_grid^2)*kron(I,DR) - spdiags(DMHZ*c_33(:),0,N_grid^2,N_grid^2)*kron(I,DL);

        %DX_12_DZ = kron(D,I)*spdiags(c_12(:),0,N_grid^2,N_grid^2)*kron(I,D);
        %DZ_33_DX = kron(I,D)*spdiags(c_33(:),0,N_grid^2,N_grid^2)*kron(D,I);

        %DZ_22_DZ = spdiags(DPHZ*c_22(:),0,N_grid^2,N_grid^2)*kron(I,DR) - spdiags(DMHZ*c_22(:),0,N_grid^2,N_grid^2)*kron(I,DL);
        %DX_33_DX = spdiags(DPHX*c_33(:),0,N_grid^2,N_grid^2)*kron(DR,I) - spdiags(DMHX*c_33(:),0,N_grid^2,N_grid^2)*kron(DL,I);

        %DZ_12_DX = kron(I,D)*spdiags(c_12(:),0,N_grid^2,N_grid^2)*kron(D,I);
        %DX_33_DZ = kron(D,I)*spdiags(c_33(:),0,N_grid^2,N_grid^2)*kron(I,D);


end

