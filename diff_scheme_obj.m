classdef  diff_scheme_obj
%diff_scheme instantiate finite difference scheme

    properties
        DX
        DZ
        DRX
        DRZ
        DLX
        DLZ
        DPHX
        DPHZ
        DMHX
        DMHZ
    end

    methods
        function diffs = diff_scheme_obj( h,N_grid )
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
            diffs.DX = kron(D,I);
            diffs.DZ = kron(I,D);
            clear D;

            % DR, DL : 1-D forward, backward f.d. ops
            e = ones(N_grid - 1,1);
            DR = (sparse(diag(e,1)) - I)/(h);
            DL = (I - sparse(diag(e,-1)))/(h);

            % DRX, DRZ, DLX, DLZ : 2-D f., b. f.d. ops
            diffs.DRX = kron(DR,I);
            diffs.DRZ = kron(I,DR);
            diffs.DLX = kron(DL,I);
            diffs.DLZ = kron(I,DL);
            clear DR;
            clear DL;

            % DPH, DMH, etc. : averaging ops
            e = ones(N_grid - 1,1);
            DPH = 0.5*(sparse(diag(e,1)) + I);
            DMH = 0.5*(sparse(diag(e,-1)) + I);
            diffs.DPHX = kron(DPH,I);
            diffs.DPHZ = kron(I,DPH);
            diffs.DMHX = kron(DMH,I);
            diffs.DMHZ = kron(I,DMH);
            clear DPH;
            clear DMH;

            clear e;
            clear I;
        end
    end

end

