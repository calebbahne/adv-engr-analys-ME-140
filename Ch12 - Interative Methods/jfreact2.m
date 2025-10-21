function [J,f]=jfreact2(x,varargin)
% jfreact2: Calculates the Jacobian matrix J and the function vector f
%           for a system of two nonlinear equations, f1(x1, x2) and f2(x1, x2).
%
% NOTE: This specific implementation uses the finite-difference approximation
%       for partial derivatives and assumes the system is defined by two
%       external functions, u(x1, x2) and v(x1, x2), where f1=u and f2=v.
%
% input:
%   x = current root estimate vector [x1; x2]
%   varargin = optional parameters to pass to u and v functions (not used here)
%
% output:
%   J = The 2x2 Jacobian matrix of partial derivatives
%   f = The 2x1 function vector [f1; f2]
%
% -------------------------------------------------------------------------
% THE JACOBIAN MATRIX (J) STRUCTURE:
%
% J = [ df1/dx1, df1/dx2 ]   -> The first row is the gradient of f1
%     [ df2/dx1, df2/dx2 ]   -> The second row is the gradient of f2
% -------------------------------------------------------------------------

% SETUP OF EXTERNAL FUNCTIONS u AND v
% -------------------------------------------------------------------------
% This function relies on the existence of two external functions, u and v,
% which represent the two nonlinear equations in the system:
% f1(x1, x2) = u(x1, x2)
% f2(x1, x2) = v(x1, x2)
%
% These functions **must** be defined separately as M-files (e.g., u.m and v.m)
% in the MATLAB search path, with the syntax:
%
% function result = u(x1, x2)
%   result = <equation for f1>;
% end
%
% **Can they be anonymous functions (@(x))?**
% While anonymous functions can be used, they are often less practical here
% because the main `jfreact2` function calls them directly by name (`u` and `v`),
% assuming they are independent function files. If using anonymous functions,
% they would need to be defined as global or passed in as additional
% arguments (varargin) to `jfreact2`, and the current code structure would
% need to be modified to handle those function handles. The standard approach
% is to use separate M-files.
% -------------------------------------------------------------------------

del=0.000001; % Define the small perturbation size (delta) for finite differencing

% --- 1. Finite-Difference Approximation for Partial Derivatives ---
% This block approximates the partial derivatives (elements of the Jacobian).
% The formula used is the forward difference: df/dx_i = (f(x_i + dx_i) - f(x_i)) / dx_i
% Note: The step size is defined relative to the current variable value: dx_i = del * x(i)

% Row 1: Partial derivatives of the first function (u or f1)
df1dx1=(u(x(1)+del*x(1),x(2))-u(x(1),x(2)))/(del*x(1)); % df1/dx1
df1dx2=(u(x(1),x(2)+del*x(2))-u(x(1),x(2)))/(del*x(2)); % df1/dx2

% Row 2: Partial derivatives of the second function (v or f2)
df2dx1=(v(x(1)+del*x(1),x(2))-v(x(1),x(2)))/(del*x(1)); % df2/dx1
df2dx2=(v(x(1),x(2)+del*x(2))-v(x(1),x(2)))/(del*x(2)); % df2/dx2

% --- 2. Construct the Jacobian Matrix ---
J=[df1dx1 df1dx2;df2dx1 df2dx2];

% --- 3. Construct the Function Vector ---
% This evaluates the functions themselves at the current estimate x.
f1=u(x(1),x(2));
f2=v(x(1),x(2));
f=[f1;f2];
end
