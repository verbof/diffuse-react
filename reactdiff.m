% REACTDIFF script for testing stencil-based CG for Diffusion/Reaction problems.
% 

%% SET PARAMS

nx      = 64;
ny      = 64;
assert(nx == ny, 'Grid must be square. Set nx = ny.');
N       = nx * ny;
x0      = 0;
xmax    = 1;
y0      = 0;
ymax    = 1;
tlim    = 0.01;
tsteps  = 10;

maxitCG = 100;
tolCG   = 1e-4;

maxitNT = 20;
tolNT   = 1e-8;

T  = linspace(0, tlim, tsteps);
dt = tlim/tsteps;

bndN = 0.00;
bndS = 0.00;
bndW = 0.00;
bndE = 0.00;

bnd.bndN = bndN;
bnd.bndS = bndS;
bnd.bndE = bndE;
bnd.bndW = bndW;

dx = (xmax-x0)/(nx);
dy = (ymax-y0)/(ny);

%D = 1;          % Diffusion constant
%R = 1;
F = 0.5;        % Initial condition

alpha = dx^2/dt;
dxs   = 100 / dx^2;


global A;
A = -gallery('poisson', nx);

%% SET INITIAL CONDITIONS
X = zeros(nx,ny);
    
% Center and radius of the circle
c  = [xmax-1/4, y0+1/4];
r  = min(c) / 2;

for j = 1:nx
    x = x0 + j*dx;
    for k = 1:ny
        p = [x, y0+k*dy];
        %disp(p);
        if norm(c-p) <= r
            X(j,k) = F;
        end
    end
end

%% SET BOUNDARY CONDITIONS

X(1,:)   = bndN;
X(end,:) = bndS;
X(:,1)   = bndW;
X(:,end) = bndE;

X = X(:);

%% RUN TIME-S

figure('Name',sprintf('%d x %d grid', nx, ny));

global f;

tTot = tic;

ts = 1;
for t = T
    
    fprintf('T-step %4d/%4d:\n', ts, tsteps);
    
    Xold = X;
    imagesc(reshape(X, nx,nx)); colorbar; pause(0.5);

    convNT = false;
    for k = 1:maxitNT
        
        %f = @(X) drstencil(X, Xold, nx, ny, alpha, dxs);
        f = @(X) diffusereact(X, Xold, alpha, dxs);
        
        X = f(Xold);
        
        fB = norm(X);
        fprintf('NWT: it=%3d, ||f(B)|| = %.4e   |  ', k, fB);
        if (fB < tolNT)
            fprintf('Newton method has converged!\n');
            convNT = true;
            break;
        end
    
        [dX, flg, relres, its] = stencil_pcg(f, X, tolCG, maxitCG);
    
        
        fprintf('CG: ');
        switch flg
        case 1
            fprintf('CG iterated %d times but did not converge.\n', maxitCG);
        case 2
            fprintf('CG preconditioner M was ill-conditioned.\n');
        case 3
            fprintf('CG stagnated.\n');
        case 4
            fprintf('CG: one of the scalar quantities is too small or too large.\n');
        end
        
        if flg > 0
            break;
        end
    
        fprintf('residual %.4e after %3d iterations.\n', relres, its);
        X    = Xold - dX;
        Xold = X;
        %spy(reshape(dX, nx,nx))
        
    end
    
    imagesc(reshape(X, nx,nx)); colorbar;
    
    if convNT
        ts = ts+1;
        continue;
    else
        error('Newton method FAILED to converge.');
    end
    
        
    
end

fprintf('-----------------------------------\n');
fprintf('Simulation took %f seconds.\n', toc(tTot));



