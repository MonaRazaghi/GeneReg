function indices = ParseCmplx(rule)
paran = strsplit(rule,'OR ');
indices = cell(size(paran,2),1);

global GeneNames;

for i = 1:size(paran,2)
    t = strrep(paran{1,i},'  ',' ');
    t = strtrim(t);
    y = find(t==' '|t=='	');
    l = length(y);
    if l >= 1
        first = t(1 : y(1) - 1);
        first = strtrim(first);
        if ~strcmpi(first,')')
            if ~strcmpi(first,'(')
                k = length(first);
                while strcmp(first(k),')')
                    first = first(1:k-1);
                    k = k - 1;
                end
                
                k = 1;
                while strcmp(first(k),'(')
                    if k + 1 <= length(first)
                        first = first(k+1:end);
                    else
                        first = '';
                        break
                    end
                    
                end
                if ~isempty(first)
                    indices{i,1} = union(indices{i,1},find(strcmp(first,GeneNames)));
                end
            end
        end
        
        for j = 1:l-1
            text = t(y(j) + 1 : y(j+1) - 1);
            text = strtrim(text);
            if ~strcmpi(text,'AND ')
                if ~strcmpi(text,')')
                    if ~strcmpi(text,'(')
                        k = length(text);
                        while strcmp(text(k),')')
                            text = text(1:k-1);
                            k = k - 1;
                        end
                        
                        k = 1;
                        while strcmp(text(k),'(')
                            text = text(k+1:end);
                        end
                       
                        indices{i,1} = union(indices{i,1},find(strcmp(text,GeneNames)));
                    end
                end
                
            end
        end
        last = t(y(l) + 1 : length(t));
        last = strtrim(last);
        if ~strcmpi(last,')')
            if ~strcmpi(last,'(')
                k = length(last);
                while strcmp(last(k),')')
                    if k - 1 >= 1 
                        last = last(1:k-1);
                        k = k - 1;
                    else
                        last = '';
                        break
                    end
                end
                
                k = 1;
                if ~isempty(last)
                    while strcmp(last(k),'(')
                        last = last(k+1:end);
                    end
                end
                
                if ~isempty(last)
                    indices{i,1} = union(indices{i,1},find(strcmp(last,GeneNames)));
                end
            end
        end
    else
        text = paran{1,i};
        text = strtrim(text);
        k = length(text);
        while strcmp(text(k),')')
            text = text(1:k-1);
            k = k - 1;
        end
        
        k = 1;
        while strcmp(text(k),'(')
            text = text(k+1:end);
        end
        
        indices{i,1} = union(indices{i,1},find(strcmp(text,GeneNames)));
    end
end
end

