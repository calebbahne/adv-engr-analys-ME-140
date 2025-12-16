function [y_eq_str, system_type, evals_str, evec1_str, evec2_str] = solve_ode_system(A)
% ANALYZE_ODE_SYSTEM_VARIABLES Solves a 2x2 linear system x' = Ax.
%   This function calculates the eigenvalues and eigenvectors, classifies the 
%   system type (e.g., Undamped, Overdamped), and generates the final 
%   symbolic displacement equation y(t) = x1(t).
%   function [y_eq, type, evals, v1, v2] = solve_ode_system(A)
%
% Inputs:
%   A: The 2x2 coefficient matrix [0 1; -k/m -c/m]
%
% Outputs:
%   system_type: Classification string (e.g., 'Pure Oscillation (Undamped)')
%   evals_str:   Formatted string of the two eigenvalues (L1 and L2)
%   evec1_str:   Formatted string of the first eigenvector (v1)
%   evec2_str:   Formatted string of the second eigenvector (v2)
%   y_eq_str:    Symbolic equation for the displacement y(t)

    % 1. INITIAL SETUP
    syms t C1 C2 % Define symbolic variables for time and integration constants
    
    % 2. EIGENVALUE AND EIGENVECTOR CALCULATION
    [V, D] = eig(A);
    lambda1 = D(1,1);
    lambda2 = D(2,2);
    v1 = V(:,1);
    v2 = V(:,2);
    
    % Extract components for easy readability and comparison
    alpha = real(lambda1);
    beta = imag(lambda1);
    
    % Initialize output variables
    system_type = 'Error: Eigenvalue case not supported (e.g., Repeated Real)';
    y_sol = 0; % Symbolic solution variable
    
    % 3. CONDITIONAL LOGIC (SYSTEM CLASSIFICATION & SOLUTION GENERATION)
    
    % --- CASE 1: Purely Imaginary (Undamped Oscillation) ---
    % Condition: alpha is near zero and beta is non-zero (e.g., +/- 2i)
    if abs(alpha) < 1e-9 && abs(beta) > 1e-9 
        system_type = 'Pure Oscillation (Undamped)';
        
        a = real(v1);
        b = imag(v1);
        
        % General formula y(t) for alpha=0: Linear combination of Real and Imaginary parts
        y_sol = C1*(a(1)*cos(beta*t) - b(1)*sin(beta*t)) + ...
                C2*(a(1)*sin(beta*t) + b(1)*cos(beta*t));

    % --- CASE 2: Complex Conjugate (Underdamped) ---
    % Condition: alpha is non-zero and beta is non-zero (e.g., -1 +/- 2i)
    elseif abs(alpha) > 1e-9 && abs(beta) > 1e-9 
        system_type = 'Damped Oscillation (Underdamped)';

        a = real(v1);
        b = imag(v1);
        
        % General formula y(t): e^(at) * [C1*(...)+C2*(...)]
        y_sol = exp(alpha*t) * ( ...
            C1 * (a(1)*cos(beta*t) - b(1)*sin(beta*t)) + ...
            C2 * (a(1)*sin(beta*t) + b(1)*cos(beta*t)) ...
        );
        
    % --- CASE 3: Real and Distinct (Overdamped) ---
    % Condition: lambda is real and L1 != L2 (e.g., -1, -4)
    elseif isreal(lambda1) && abs(lambda1 - lambda2) > 1e-9
        system_type = 'Non-oscillatory Decay (Overdamped)';

        % General formula y(t): C1*v1(1)*exp(l1*t) + C2*v2(1)*exp(l2*t)
        y_sol = C1*v1(1)*exp(lambda1*t) + C2*v2(1)*exp(lambda2*t);

    end

    % 4. FORMATTING OUTPUT STRINGS
    
    % Use num2str and sprintf for clean, readable complex number display
    % Eigenvalues
    evals_str = sprintf('L1: %s \nL2: %s', num2str(lambda1, 5), num2str(lambda2, 5));
    
    % Eigenvector 1: [v1(1) ; v1(2)]
    evec1_str = sprintf('[%s ; %s]', num2str(v1(1), 5), num2str(v1(2), 5));
    
    % Eigenvector 2: [v2(1) ; v2(2)]
    evec2_str = sprintf('[%s ; %s]', num2str(v2(1), 5), num2str(v2(2), 5));
    
    % Final Y Equation
    y_eq_str = char(simplify(y_sol));
    
end