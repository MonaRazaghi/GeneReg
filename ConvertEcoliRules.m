function ConvertedRules = ConvertEcoliRules(rules,gene_names)

ruleNum = size(rules,1);

for i = 1:ruleNum

    if ~isempty(rules{i})
        t = char(rules{i});
        x_pos = strfind(t,'x');
        while ~isempty(x_pos)
            k = x_pos(1);
            while t(k) ~= ')'
                k = k + 1;
            end
            t = strrep(t,t(x_pos(1):k),gene_names{str2num(t(x_pos(1)+2:k-1))});
            x_pos = [];
            x_pos = strfind(t,'x');
        end
        t = strrep(t,'&','and');
        t = strrep(t,'|','or');
        
        ConvertedRules{i,1} = char(t);
    end
end


end