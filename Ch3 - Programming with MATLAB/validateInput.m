function validateInput(input, expectedType, range)
%VALIDATEINPUT: checks type (& optionally the range) of an input
% Throws an error if invalid.
%
% validateInput(input, expectedType)
% validateInput(input, expectedType, range)
% inputs:
%   input           = variable to validate
%   expectedType    = 'numeric', 'char', 'logical'
%   range           = optional [MIN MAX], only if 'numeric'

if ~isa(input, expectedType)
    error('Input must be of type %s', expectedType);
end

if nargin == 3 && ~isempty(range)
    if input < range(1) || input > range(2)
        error('Input must be between %g and %g', range(1), range(2));
    end
end