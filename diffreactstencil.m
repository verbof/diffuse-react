function S = diffreactstencil(U, x_old, nx, ny, bnd, alpha, dxs)
% DIFFREACTSTENCIL computes stencil as inner + bounds + corners

bndN = bnd.bndN;
bndS = bnd.bndS;
bndE = bnd.bndE;
bndW = bnd.bndW;

S = zeros(size(U));

%% Interior grid points
for j = 2:ny-1          %(int j=1; j < jend; j++) {
    for i = 2:nx-1      %(int i=1; i < iend; i++) {
        S(i,j) = -(4.0 + alpha) * U(i,j) ...               % central point
                                + U(i-1,j) + U(i+1,j) ... % east and west
                                + U(i,j-1) + U(i,j+1) ... % north and south
                                + alpha * x_old(i,j)  ... 
                                + dxs * U(i,j) * (1.0 - U(i,j));
    end
end
% end Interior grid points

%% East boundary
i = nx;                 %int i = nx - 1;
for j = 2:ny-1          %(int j = 1; j < jend; j++)
    S(i,j) = -(4.0 + alpha) * U(i,j)     ...
                + U(i-1,j) + U(i,j-1) + U(i,j+1) ...
                + alpha*x_old(i,j) + bndE(j)     ...
                + dxs * U(i,j) * (1.0 - U(i,j));
end
% end East boundary

%% West boundary
i = 1;                  %int i = 0;
for j = 2:ny-1 %(int j = 1; j < jend; j++)
    S(i,j) = -(4.0 + alpha) * U(i,j)         ...
                + U(i+1,j) + U(i,j-1) + U(i,j+1) ...
                + alpha * x_old(i,j) + bndW(j)   ...
                + dxs * U(i,j) * (1.0 - U(i,j));
end
% end West boundary

%% North boundary, plus NE and NW corners
j = ny;           %int j = ny - 1;

% NW corner
i = 1;            %int i = 0; // NW corner
S(i,j) = -(4.0 + alpha) * U(i,j)         ...
            + U(i+1,j) + U(i,j-1)       ...
            + alpha * x_old(i,j) + bndW(j) + bndN(i) ...
            + dxs * U(i,j) * (1.0 - U(i,j));
% end NW corner

% North boundary w/o corners
for i = 2:nx-1        %(int i = 1; i < iend; i++)
    S(i,j) = -(4.0 + alpha) * U(i,j)        ...
                + U(i-1,j) + U(i+1,j) + U(i,j-1) ...
                + alpha*x_old(i,j) + bndN(i)    ...
                + dxs * U(i,j) * (1.0 - U(i,j));
end
% end North boundary w/o corners

% NE corner
i = nx;        % int i = nx-1; // NE corner
S(i,j) = -(4.0 + alpha) * U(i,j)        ...
            + U(i-1,j) + U(i,j-1)       ...
            + alpha * x_old(i,j) + bndE(j) + bndN(i)    ...
            + dxs * U(i,j) * (1.0 - U(i,j));
% end NE corner

%% South boundary, plus SW and SE corners
j = 1;        % int j = 0;
   
% SW corner
i = 1;        % int i = 0; // SW corner
S(i,j) = -(4.0 + alpha) * U(i,j)    ...
            + U(i+1,j) + U(i,j+1)   ...
            + alpha * x_old(i,j) + bndW(j) + bndS(i) ...
            + dxs * U(i,j) * (1.0 - U(i,j));
% end SW corner

% South boundary w/o corners
for i = 2:nx-1      %(int i = 1; i < iend; i++)
    S(i,j) = -(4.0 + alpha) * U(i,j)                ...
            + U(i-1,j) + U(i+1,j) + U(i,j+1)    ...
            + alpha * x_old(i,j) + bndS(i)      ...
            + dxs * U(i,j) * (1.0 - U(i,j));
end
% end South boundary w/o corners

% SE corner
i = nx;     % int i = nx - 1; // SE corner
S(i,j) = -(4.0 + alpha) * U(i,j)        ...
            + U(i-1,j) + U(i,j+1)       ...
            + alpha * x_old(i,j) + bndE(j) + bndS(i)    ...
            + dxs * U(i,j) * (1.0 - U(i,j));
% end SW corner

end




