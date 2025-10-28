function [eigval, eigvect] = powereig(A, es, maxit)
% powereig - Computes the dominant eigenvalue and eigenvector using the Power Method
%
% Syntax:
%   [eigval, eigvect] = powereig(A, es, maxit)
%
% Inputs:
%   A     - Square matrix for which the dominant eigenvalue/vector is sought
%   es    - Desired relative error tolerance (stopping criterion, e.g. 0.0001)
%   maxit - Maximum number of iterations to perform
%
% Outputs:
%   eigval  - Dominant (largest magnitude) eigenvalue of A
%   eigvect - Corresponding normalized eigenvector
%
% Example:
%   A = [4 1; 2 3];
%   [eigval, eigvect] = powereig(A, 0.0001, 100)
%
% Description:
%   This function implements the Power Method, an iterative algorithm used 
%   to find the largest eigenvalue (in magnitude) and its associated eigenvector.
%   The algorithm repeatedly multiplies an arbitrary vector by the matrix A,
%   normalizes the result, and estimates convergence based on the change in eigenvalue.

    n = length(A);              % Determine the size of matrix A
    eigvect = ones(n,1);          % Initialize eigenvector guess (all ones)
    eigval = 1;                   % Initial guess for eigenvalue
    iter = 0;                   % Iteration counter
    ea = 100;                   % Initialize approximate relative error (%)

    while (1)
        eigvalold = eigval;         % Save previous eigenvalue estimate
        
        eigvect = A * eigvect;      % Multiply A by current eigenvector estimate
        eigval = max(abs(eigvect)); % Approximate new eigenvalue (dominant component)
        eigvect = eigvect ./ eigval;  % Normalize eigenvector
        
        iter = iter + 1;        % Increment iteration counter
        
        if eigval ~= 0
            ea = abs((eigval - eigvalold) / eigval) * 100; % Compute approximate relative error
        end
        
        % Stop if error is below threshold or max iterations reached
        if ea <= es || iter >= maxit
            break
        end
    end
end
