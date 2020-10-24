clear;
clc;

load iJR904
model = iJR904;

[modelIrrev,~,~,irrev2rev] = convertToIrreversible(model);
[mNum,rNum] = size(modelIrrev.S);

biomass = find(modelIrrev.c == 1); % 
Obj = 329; % Ethanol exchange reaction
biomass_min = 0.47; % 50% of maximum biomass
Mu_max = 1000;
modelIrrev.c = zeros(rNum,1);

GeneRules = modelIrrev.grRules;

GeneNum = length(modelIrrev.genes);
max_deletion = GeneNum;

global GeneNames
GeneNames = modelIrrev.genes;

Enzymes = [];
i_E = 0;
for r = 1:length(GeneRules)
    if ~isempty(GeneRules{r})
        r_enzymes = FindEnzymes(GeneRules{r});
        for i = 1:size(r_enzymes,1)
            found = 0;
            for j = 1:size(Enzymes,1)
                R = r_enzymes{i,:};
                E = Enzymes{j,:}; 
                if length(R) == length(E)
                    if sum(R == E) == length(R)
                        found = 1;
                        break
                    end
                end
            end
            if found == 0
                i_E = i_E + 1;
                Enzymes{i_E,1} = r_enzymes{i,:};
            end
        end
    end
end

Gene_Enzyme = zeros(length(GeneNames),length(Enzymes));
for i = 1:length(Enzymes)
    E = Enzymes{i,:};
    Gene_Enzyme(E,i) = 1;
end

Enzyme_Reaction = zeros(length(Enzymes),rNum);
for r = 1:rNum
    if ~isempty(GeneRules{r})
        r_enzymes = FindEnzymes(GeneRules{r});
        for i = 1:size(r_enzymes,1)
            found = 0;
            for j = 1:size(Enzymes,1)
                R = r_enzymes{i,:};
                E = Enzymes{j,:}; 
                if length(R) == length(E)
                    if sum(R == E) == length(R)
                        found = j;
                        break
                    end
                end
            end
            if found > 0
                Enzyme_Reaction(j,r) = 1;
            end
        end
    end
end

eNum = length(Enzymes);
%%
inx_inequality = 1;

% 4
A_inequality(inx_inequality,2*rNum+eNum+1:2*rNum+eNum+GeneNum) = -1;
b_inequality(inx_inequality,1) = max_deletion - GeneNum;

% 5
inx_inequality = inx_inequality + 1;
A_inequality(inx_inequality,biomass) = -1;
b_inequality(inx_inequality,1) = -biomass_min;

% 6
for r = 1:rNum
    inx_inequality = inx_inequality + 1;
    A_inequality(inx_inequality,2*rNum+eNum+GeneNum+r) = 1;
    A_inequality(inx_inequality,rNum+r) = Mu_max;
    b_inequality(inx_inequality,1) = Mu_max;
end

% 7.1
for r = 1:rNum
    inx_inequality = inx_inequality + 1;
    A_inequality(inx_inequality,r) = -1;
    A_inequality(inx_inequality,rNum+r) = modelIrrev.lb(r);
    b_inequality(inx_inequality,1) = 0;
end

% 7.2
for r = 1:rNum
    inx_inequality = inx_inequality + 1;
    A_inequality(inx_inequality,r) = 1;
    A_inequality(inx_inequality,rNum+r) = -modelIrrev.ub(r);
    b_inequality(inx_inequality,1) = 0;
end

% 8
for r = 1:rNum
    if sum(Enzyme_Reaction(:,r)) > 0
        e_inx = find(Enzyme_Reaction(:,r) == 1);
        for i = 1:length(e_inx)
            inx_inequality = inx_inequality + 1;
            A_inequality(inx_inequality,rNum+r) = -1;
            A_inequality(inx_inequality,2*rNum+e_inx(i)) = 1;
            b_inequality(inx_inequality,1) = 0;
        end
    end
end


% 9
for r = 1:rNum
    if sum(Enzyme_Reaction(:,r)) > 0
        inx_inequality = inx_inequality + 1;
        A_inequality(inx_inequality,rNum+r) = 1;
        A_inequality(inx_inequality,2*rNum+find(Enzyme_Reaction(:,r) == 1)) = -1;
        b_inequality(inx_inequality,1) = 0;
    end
end

% 10
for e = 1:eNum
    if sum(Gene_Enzyme(:,e)) > 0
        inx_inequality = inx_inequality + 1;
        A_inequality(inx_inequality,2*rNum+e) = -1;
        A_inequality(inx_inequality,2*rNum+eNum+find(Gene_Enzyme(:,e) == 1)) = 1;
        b_inequality(inx_inequality,1) = sum(Gene_Enzyme(:,e)) - 1;
    end
end

% 11
for e = 1:eNum
    if sum(Gene_Enzyme(:,e)) > 0
        g_inx = find(Gene_Enzyme(:,e) == 1);
        for i = 1:length(g_inx)
            inx_inequality = inx_inequality + 1;
            A_inequality(inx_inequality,2*rNum+e) = 1;
            A_inequality(inx_inequality,2*rNum+eNum+g_inx(i)) = -1;
            b_inequality(inx_inequality,1) = 0;
        end
    end
end

%%
inx_equality = 1;

A_equality(inx_equality:mNum,1:rNum) = modelIrrev.S;
b_equality(inx_equality:mNum,1) = zeros(mNum,1);

inx_equality = mNum;

% 2
inx_equality = inx_equality + 1;
A_equality(inx_equality,3*rNum+eNum+GeneNum+1:3*rNum+eNum+GeneNum+mNum) = (modelIrrev.S(:,biomass))';
A_equality(inx_equality,2*rNum+eNum+GeneNum+biomass) = 1;
b_equality(inx_equality,1) = 1;

% 3
for r = 1:rNum
    if r ~= biomass
        inx_equality = inx_equality + 1;
        A_equality(inx_equality,3*rNum+eNum+GeneNum+1:3*rNum+eNum+GeneNum+mNum) = (modelIrrev.S(:,r))';
        A_equality(inx_equality,2*rNum+eNum+GeneNum+r) = 1;
        b_equality(inx_equality,1) = 0;        
    end
end

% 12
inx_equality = inx_equality + 1;
A_equality(inx_equality,biomass) = 1;
A_equality(inx_equality,2*rNum+eNum+GeneNum+biomass) = -biomass_min;
b_equality(inx_equality,1) = 0;

%
A_inequality(:,3*rNum+eNum+GeneNum+1:3*rNum+eNum+GeneNum+mNum) = zeros(size(A_inequality,1),mNum);


%%
LB = zeros(3*rNum+eNum+GeneNum+mNum,1);
UB = ones(3*rNum+eNum+GeneNum+mNum,1);

LB(1:rNum) = modelIrrev.lb;
UB(1:rNum) = modelIrrev.ub;

LB(2*rNum+eNum+GeneNum+1:3*rNum+eNum+GeneNum+mNum) = -1000;
UB(2*rNum+eNum+GeneNum+1:3*rNum+eNum+GeneNum+mNum) = 1000;

LB(2*rNum+eNum+1:2*rNum+eNum+GeneNum) = 1;
%%
alpha = 1;
f = zeros(size(A_equality,2),1);
f(Obj) = -1;
f(2*rNum+eNum+1:2*rNum+eNum+GeneNum) = -alpha;
intcon = rNum+1:rNum+eNum+GeneNum;
options = optimoptions('intlinprog','IntegerPreprocess','none','HeuristicsMaxNodes', 100,'MaxTime',345600,'Heuristics','round','BranchRule','maxpscost','LPMaxIterations',3e5,'IntegerPreprocess', 'none','LPPreprocess','none');

[x,fval,~] = intlinprog(f,intcon,A_inequality,b_inequality,A_equality,b_equality,LB,UB,options);

x(biomass)
x(Obj)

numel(find(x(2*rNum+eNum+1:2*rNum+eNum+GeneNum) == 0))
% GeneNames(find(x(2*rNum+eNum+1:2*rNum+eNum+GeneNum) == 0))