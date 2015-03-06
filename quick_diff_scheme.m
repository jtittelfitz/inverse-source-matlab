function [DX,DZ] = quick_diff_scheme( h,N_grid )
%diff_scheme instantiate finite difference scheme

    I = speye(N_grid);

	% D : 1-D centered finite difference operator
    e = (1/2)*ones(N_grid - 1,1);
    D = sparse(diag(e,1));
    D = D - D';
    %D(1,1:2) = [-1,1];
    %D(N_grid,N_grid - 1:N_grid) = [-1,1];

    D(1,1) = -3/2;
    D(1,2) = 4/2;
    D(1,3) = -1/2;
    D(N_grid, N_grid - 2) = 1/2;
    D(N_grid, N_grid - 1) = -4/2;
    D(N_grid, N_grid) = 3/2;

    D = D/h;

	% DX, DZ : 2-D center f.d. ops
    DX = kron(D,I);
    DZ = kron(I,D);    
end

