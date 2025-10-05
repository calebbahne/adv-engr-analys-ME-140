function [x1,ea,iter] = fixpt_autog(f,x0,es,maxit)

% NOT CURRENTLY WORKING
% fixpt_autog: automatic fixed-point iteration
%   [x1,ea,iter] = fixpt_autog(f,x0,es,maxit)
%   Automatically rearranges f(x)=0 into x=g(x),
%   converts g(x) to a function handle, and performs fixed-point iteration.
%
% Input:
%   f     = function handle for f(x)
%   x0    = initial guess
%   es    = stopping criterion for approximate relative error (%)
%           (default = 1e-6 %)
%   maxit = maximum iterations (default = 50)
%
% Output:
%   x1    = estimated root
%   ea    = approximate relative error (%)
%   iter  = number of iterations performed

    % --- Input handling ---
    if nargin < 2
        error('At least f and x0 must be provided.');
    end
    if nargin < 3 || isempty(es)
        es = 1e-6;
    end
    if nargin < 4 || isempty(maxit)
        maxit = 50;
    end

    % --- Step 1: Convert f to symbolic form ---
    syms x
    fSym = f(x);              % Convert numeric function handle f into a symbolic expression
    eqn = fSym == 0;          % Create the symbolic equation f(x) = 0

    gSym = rhs(isolate(eqn, x)); % Rearrange equation to isolate x, then take the right-hand side (g(x))

    % --- Step 3: Convert g(x) back to function handle ---
    g = matlabFunction(gSym);

    % --- Step 4: Run fixed-point iteration (internal loop) ---
    iter = 0;
    ea = 100; % start with 100% error
    while true
        x1 = g(x0);
        iter = iter + 1;
        if x1 ~= 0
            ea = abs((x1 - x0) / x1) * 100;
        end
        if (ea <= es) || (iter >= maxit)
            break;
        end
        x0 = x1;
    end
end
