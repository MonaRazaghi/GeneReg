function [GeneList] = extractingGenes2(t)
s = find(t==' '|t=='	');
GeneList{1,1} = strtrim(t(1:s(1)-1));
ii = 1;
for i = 1:length(s)-1
    if ~isempty(strtrim(t(s(i)+1:s(i+1)-1)))
        if ~strcmp('and',strtrim(t(s(i)+1:s(i+1)-1)))
            if ~strcmp('or',strtrim(t(s(i)+1:s(i+1)-1)))
                ii = ii + 1;
                GeneList{ii,1} = strtrim(t(s(i)+1:s(i+1)-1));
            end
        end
    end
end

ii = ii + 1;
GeneList{ii,1} = strtrim(t(s(length(s))+1:length(t)));
end