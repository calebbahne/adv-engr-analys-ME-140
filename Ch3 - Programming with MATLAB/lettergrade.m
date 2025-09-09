function grade = lettergrade(score)
%grade Returns letter grade given numeric grade
% score = numeric grade

% check for bad inputs
validateInput(score,'numeric',[0 100]);

% Can't do switch, cuz this is range -> category
if score >= 90
    grade = 'A';
elseif score >= 80
    grade = 'B';
elseif score >= 70
    grade = 'C';
elseif score >= 60
    grade = 'D';
else
    grade = 'F';
end
end