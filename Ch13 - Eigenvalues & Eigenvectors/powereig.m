function [eval, evect] = powereig(A, es, maxit)
% powereig - Computes the dominant eigenvalue and eigenvector using the Power Method
%
% Syntax:
%   [eval, evect] = powereig(A, es, maxit)
%
% Inputs:
%   A     - Square matrix for which the dominant eigenvalue/vector is sought
%   es    - Desired relative error tolerance (stopping criterion, e.g. 0.0001)
%   maxit - Maximum number of iterations to perform
%
% Outputs:
%   eval  - Dominant (largest magnitude) eigenvalue of A
%   evect - Corresponding normalized eigenvector
%
% Example:
%   A = [4 1; 2 3];
%   [eval, evect] = powereig(A, 0.0001, 100)
%
% Description:
%   This function implements the Power Method, an iterative algorithm used 
%   to find the largest eigenvalue (in magnitude) and its associated eigenvector.
%   The algorithm repeatedly multiplies an arbitrary vector by the matrix A,
%   normalizes the result, and estimates convergence based on the change in eigenvalue.

    n = length(A);              % Determine the size of matrix A
    evect = ones(n,1);          % Initialize eigenvector guess (all ones)
    eval = 1;                   % Initial guess for eigenvalue
    iter = 0;                   % Iteration counter
    ea = 100;                   % Initialize approximate relative error (%)

    while (1)
        evalold = eval;         % Save previous eigenvalue estimate
        
        evect = A * evect;      % Multiply A by current eigenvector estimate
        eval = max(abs(evect)); % Approximate new eigenvalue (dominant component)
        evect = evect ./ eval;  % Normalize eigenvector
        
        iter = iter + 1;        % Increment iteration counter
        
        if eval ~= 0
            ea = abs((eval - evalold) / eval) * 100; % Compute approximate relative error
        end
        
        % Stop if error is below threshold or max iterations reached
        if ea <= es || iter >= maxit
            break
        end
    end
end
