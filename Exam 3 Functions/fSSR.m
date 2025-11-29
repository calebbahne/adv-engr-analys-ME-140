function f = fSSR(a,xm,ym)
%fSSR: f = fSSR(a,xm,ym)
% F = a0*v^a1
% a(1) = a0, a(2) = a1
% sum of squares

% call afterwards****
% a = fminsearch(@fSSR, [1 1], [], v, F)
%   [1 1] is an initial guess for [a0, a1], [] is a placeholder for options
    F_pred = a(1) * xm.^a(2);  % change
    % F_pred = a(1)*xm.*exp(a(2)*xm);
    f = sum((ym-F_pred).^2);  % change the function
end