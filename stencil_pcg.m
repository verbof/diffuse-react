function [x, flg, resNew, its] = stencil_pcg(f, b, tol, maxit)
% STENCIL_PCG computes CG with 1st order approximation for f

    its = 0;
    flg = 1;

    e = 1e-08;
    
    x0 = zeros(size(b));     % initial guess
    x  = x0;
    
    
    Fx  = f(x);
    xx  = (1+e)*x;
    Fxx = f(xx);
    r = b - (1/e) * (Fxx - Fx);         % A*x \approx (1/e) * (Fxx-Fx)
    p = r;
    
    resOld = r'*r;
    resNew = resOld;
    
    if (sqrt(resOld) < tol)
        flg = 0;
        %fprintf('STENCIL_CG: converged in 0 iterations (initial guess is good enough!\n');
        return;
    end
    
    for its = 1:maxit
        
        %fprintf('STENCIL_CG: %3d ', its);
       
        Fx = f(x);
        xx  = x + e*p;
        Fxx = f(xx);
        Ap  = (1/e) * (Fxx-Fx);     % A*p approx. by Ap
        %fprintf(' * RES. FOR Ap = %.4e * ', norm(Ap-A*p)/norm(A*p));
        
        alpha = r'*r / (p'*Ap);      
        
        x    = x + alpha * p;
        rnew = r - alpha * Ap;
        
        resNew = rnew' * rnew;
        if (sqrt(resNew) < tol)
            flg = 0;
            break;
        end

        beta = resNew / resOld;   % res=<r(k-1),r(k-1)>, resNew=<r(k+1),r(k+1)>
        p    = rnew + beta * p;
 
        r      = rnew;
        resOld = resNew;
        %fprintf(' - RELRES: %.4e\n', sqrt(resNew));
        
    end
    
    if flg == 0
        %fprintf('STENCIL_CG: converged in %d iterations.\n', its);
    else
        %fprintf('STENCIL_CG: could NOT converge.\n');
    end
    
end



