function B = diffusereact(X, Xold, alpha, dxs)
% DIFFUSEREACT computes Diffusion/Reaction using Laplacian matrix.
    
    global A;
    
    aI = alpha*speye(size(A));
    B  = (A - aI)*X + alpha*Xold + dxs * (X - X.^2);
end

