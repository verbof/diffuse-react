function S = drstencil(X, Xold, nx, ny, alpha, dxs)
% DRSTENCIL computes the stencil for a grid including BCs

    X    = reshape(X,nx,ny);
    Xold = reshape(Xold,nx,ny);
    S    = zeros(size(X));

    for j = 2:ny-1
        for i = 2:nx-1
            S(i,j) = -(4.0 + alpha) * X(i,j) ...               % central point
                                    + X(i-1,j) + X(i+1,j) ...  % east and west
                                    + X(i,j-1) + X(i,j+1) ...  % north and south
                                    + alpha * Xold(i,j)  ... 
                                    + dxs * X(i,j) * (1.0 - X(i,j));
        end
    end

    S = S(:);

end


