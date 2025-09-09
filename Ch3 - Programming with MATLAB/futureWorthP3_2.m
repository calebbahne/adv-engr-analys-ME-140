function futureWorthP3_2(P, i, n)
%futureworth: yield after compound interest, in table form
% Input
%   P = initial investment
%   i = interest rate
%   n = number of periods
% Output
%   Future worth of investments from year 1 to n
%   Table with headings for year and worth

nn = 0:n;
F = P*(1+i).^nn; % do it for every year
y = [nn;F]; % horiz table
fprintf('\n Year \tFuture Worth ($)\n');
fprintf('%5d %14.f\n', y);

end