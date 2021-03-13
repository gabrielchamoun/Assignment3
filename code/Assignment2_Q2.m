function [V] = Assignment2_Q2(nx, ny, boxL, boxW, sigma)
    
    % Conductivity Map
    global cMap
    cMap = ones(nx,ny);
    cMap(nx/2 - boxL/2:nx/2 + boxL/2,1:boxW) = sigma;
    cMap(nx/2 - boxL/2:nx/2 + boxL/2,ny-boxW:ny) = sigma;
    
    G = sparse(nx*ny,ny*nx);
    F = zeros(nx*ny,1);
    for i = 1:nx
        for j = 1:ny

            % Node Mapping
            n = j + (i-1) * ny;     % middle
            nxm = j + (i-2) * ny;	% right
            nxp = j + i * ny;       % left
            nym = j-1 + (i-1) * ny; % top 
            nyp = j+1 + (i-1) * ny; % down

            if i == 1
                G(n,n) = 1;
                F(n,1) = 1;
            elseif i == nx
                G(n,n) = 1;
            elseif j == 1
                rxm = (cMap(i,j) + cMap(i-1,j)) / 2;
                rxp = (cMap(i,j) + cMap(i+1,j)) / 2;
                ryp = (cMap(i,j) + cMap(i,j+1)) / 2;

                G(n,n) = -(rxm+rxp+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nyp) = ryp;
            elseif j == ny
                rxm = (cMap(i,j) + cMap(i-1,j)) / 2;
                rxp = (cMap(i,j) + cMap(i+1,j)) / 2;
                rym = (cMap(i,j) + cMap(i,j-1)) / 2;

                G(n,n) = -(rxm+rxp+rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
            else
                rxm = (cMap(i,j) + cMap(i-1,j)) / 2;
                rxp = (cMap(i,j) + cMap(i+1,j)) / 2;
                rym = (cMap(i,j) + cMap(i,j-1)) / 2;
                ryp = (cMap(i,j) + cMap(i,j+1)) / 2;

                G(n,n) = -(rxm+rxp+rym+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
                G(n,nyp) = ryp;
            end

        end
    end

    % Solving the set of linear equations to obtain voltage values
    dA = decomposition(G,'lu');
    V = dA\F;   % Return vector of potentials
end