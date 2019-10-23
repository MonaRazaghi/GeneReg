function Check = checkRules(rules)
N = length(rules);
Check = zeros(N,1);

for i = 1:N
    if isempty(rules{i})
        Check(i) = 0;
    else
        t = strtrim(char(rules{i}));
        AND = strfind(t,'and');
        OR = strfind(t,'or');
        if length(AND) + length(OR) == 0
            Check(i) = 1;
        elseif length(AND) == 0 
            Check(i) = 2;
        elseif length(OR) == 0 
            Check(i) = 3;
        else
            Check(i) = 4;
        end           
        
    end
end