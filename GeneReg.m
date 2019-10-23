
tic
clear;
clc;

%reading the model
model = readCbModel

%split reversible reactions
[modelIrrev,~,~,irrev2rev] = convertToIrreversible(model);

[mNum,rNum] = size(modelIrrev.S);

biomass = 13; % biomass reaction
Obj = 39;     % Succinate exchange

per_biomass = 0.2;
f_Obj = 0.3;
CC = 0.01;
L = 4*rNum;

% v_min the minimum possible flux value
% v_max the maximum possible flux value
% v0_l the minimum possible flux value in the wild type
% v0_u the maximum possible flux value in the wild type
% in case of big models it is recommended to cacaulate v0_l,v0_u, v_min and
% v_max in advance

changeCobraSolver('matlab');
OPT = optimizeCbModel(modelIrrev);

modelbiomass = changeRxnBounds(modelIrrev,char(modelIrrev.rxns{find(modelIrrev.c ~= 0)}),OPT.f,'b');
modelbiomass.c(find(modelbiomass.c ~= 0)) = 0;
[v0_l, v0_u] = fluxVariability(modelbiomass);

modelIrrev.c(find(modelIrrev.c ~= 0)) = 0;
[v_min, v_max] = fluxVariability(modelIrrev);

%% finding reversible reactions and their matches
info_rev = zeros(rNum,2); %info_rev(:,1) = 0, 1, -1: irrev, f, b    info_rev(:,2) = the matching reaction for reversible ones

for i = 1:rNum
    index = find(irrev2rev == irrev2rev(i));
    if length(index) == 2
        info_rev(index(1),1) = 1;
        info_rev(index(2),1) = -1;
        
        info_rev(index(1),2) = index(2);
        info_rev(index(2),2) = index(1);
    end
end
num_rev = nnz(info_rev(:,1))/2;

%%
% in models where there is no GPR rules in terms of genes but in terms of
% x(i) this function will be used
% modelIrrev.grRules = ConvertEcoliRules(modelIrrev.rules,modelIrrev.genes);

GeneRules = modelIrrev.grRules;
ChkList = Rules_type(GeneRules); %simple vs complex reactions in term of gene rules 0: no rule, 1: single gene, 2: isoenzyme, 3: complex enzyme, 4: complex isoenzyme

%%
GeneNum = length(modelIrrev.genes);

global GeneNames
GeneNames = modelIrrev.genes;

for i = 1:rNum
    if v_max(i) < v0_u(i)
        v0_u(i) = v_max(i);
    end
end

dd = ((1 - CC) * v0_l) + (CC * v_min);
uu = ((1 - CC) * v0_u) + (CC * v_max);

% for numerical problems, we make the following changes
v_max = v_max + 10^(-3);
v_min = v_min - 10^(-3);

%% inequality matrix

Aineq = sparse((5*rNum + 2*num_rev + 1 + GeneNum), 3 * rNum + 2 * GeneNum);
bineq = zeros((5*rNum + 2*num_rev + 1 + GeneNum), 1);

% inequalities 5-8
Aineq(1:rNum,1:rNum) = eye(rNum);
Aineq(1:rNum,rNum+1:2*rNum) = diag(dd-v_max);
bineq(1:rNum,1) = dd;

Aineq(rNum+1:2*rNum,1:rNum) = -eye(rNum);
Aineq(rNum+1:2*rNum,rNum+1:2*rNum) = diag(v0_l-v_min);
bineq(rNum+1:2*rNum,1) = -v_min;

Aineq(2*rNum+1:3*rNum,1:rNum) = eye(rNum);
Aineq(2*rNum+1:3*rNum,2*rNum+1:3*rNum) = diag(v_max-v0_u);
bineq(2*rNum+1:3*rNum,1) = v_max;

Aineq(3*rNum+1:4*rNum,1:rNum) = -eye(rNum);
Aineq(3*rNum+1:4*rNum,2*rNum+1:3*rNum) = diag(v_min-uu);
bineq(3*rNum+1:4*rNum,1) = -uu;

% inequality 9
Aineq(4*rNum+1:5*rNum,rNum+1:2*rNum) = -eye(rNum);
Aineq(4*rNum+1:5*rNum,2*rNum+1:3*rNum) = -eye(rNum);
bineq(4*rNum+1:5*rNum,1) = -(diag(eye(rNum)));

% inequality 10
Aineq(5*rNum+1,rNum+1:3*rNum) = -1;
bineq(5*rNum+1,1) = L - 2*rNum;

% inequalities 11 and 12
Atemp = zeros(num_rev,rNum);
inx = find(info_rev(:,1) == 1);
for i = 1:length(inx)
    Atemp(i,inx(i)) = -1;
    Atemp(i,info_rev(inx(i),2)) = -1;
end

Aineq(5*rNum+2:5*rNum+num_rev+1,rNum+1:2*rNum) = Atemp;
bineq(5*rNum+2:5*rNum+num_rev+1,1) = -diag(eye(num_rev));

Aineq(5*rNum+num_rev+2:5*rNum+2*num_rev+1,2*rNum+1:3*rNum) = Atemp;
bineq(5*rNum+num_rev+2:5*rNum+2*num_rev+1,1) = -diag(eye(num_rev));

% inequality 23
Aineq(1:size(Aineq,1),3*rNum+1:3*rNum+2*GeneNum) = 0;
%
for g = 1: GeneNum
    Aineq(5*rNum+2*num_rev+1+g,3*rNum+g) = 1;
    Aineq(5*rNum+2*num_rev+1+g,3*rNum+GeneNum+g) = 1;
    bineq(5*rNum+2*num_rev+1+g,1) = 1;
end

%% LB and UB

LB(1:rNum,1) = modelIrrev.lb;
LB(rNum+1:3*rNum+2*GeneNum) = zeros(2*rNum+2*GeneNum,1);

LB(biomass) = per_biomass * v0_u(biomass); % inequality 3
LB(Obj) = f_Obj * v_max(Obj); %  inequality 4

UB(1:rNum,1) = modelIrrev.ub;

% binary variables
UB(rNum+1:3*rNum+2*GeneNum,1) = ones(2*rNum+2*GeneNum,1);

%% relating genes to reactions by linear constraints

r_Aineq = size(Aineq,1); % last row index
c_Aineq = size(Aineq,2); % last col index
inx_LB = size(LB,1);

max_complexes = 100;
variable = zeros(GeneNum,max_complexes);
variable_end = 0;

for r = 1:length(GeneRules)
    
    % isoenzymes: inequalities 15 and 16
    if ChkList(r) == 2
        t = erase(char(GeneRules{r}),"(");
        t = erase(t,")");
        t = strtrim(t);
        
        genes = extractingGenes(t);
        
        %
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = 1;
        end
        Aineq(r_Aineq,rNum + r) = length(genes);
        bineq(r_Aineq,1) = 2*length(genes)-1;
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = -1;
        end
        Aineq(r_Aineq,rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1 * length(genes);
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = 1;
        end
        Aineq(r_Aineq,2*rNum + r) = length(genes);
        bineq(r_Aineq,1) = length(genes);
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = -1;
        end
        Aineq(r_Aineq,2*rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1;
        
    end
    
    % complex enzymes: inequalities 17 and 18
    if ChkList(r) == 3
        t = erase(char(GeneRules{r}),"(");
        t = erase(t,")");
        t = strtrim(t);
        genes = extractingGenes(t);
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = 1;
        end
        Aineq(r_Aineq,rNum + r) = length(genes);
        bineq(r_Aineq,1) = length(genes);
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = -1;
        end
        Aineq(r_Aineq,rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1;
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = 1;
        end
        Aineq(r_Aineq,2*rNum + r) = length(genes);
        bineq(r_Aineq,1) = 2*length(genes)-1;
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = -1;
        end
        Aineq(r_Aineq,2*rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1 * length(genes);
        
    end
    
    % down-regulation of complex isoenzymes
    if ChkList(r) == 4
        complexes = ParseCmplx(char(GeneRules{r}));
        current_G = [];
        for i = 1:size(complexes,1)
            genes = complexes{i,1};
            cmplx_inx = 0;
            for j = 1:variable_end
                if sum(variable(genes,j)) == length(genes)
                    if sum(variable(setdiff(1:GeneNum,genes),j)) == 0
                        cmplx_inx = j;
                        break
                    end
                end
            end
            if cmplx_inx == 0
                variable_end = variable_end + 1;
                variable(genes,variable_end) = 1;
                
                cmplx_inx = variable_end;
            end
            current_G = union(current_G,cmplx_inx);
            
            % 
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + genes(j)) = 1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*cmplx_inx - 1)) = -1 * length(genes);
            bineq(r_Aineq,1) = 0;
            
            % 
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + genes(j)) = -1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*cmplx_inx - 1)) = length(genes);
            bineq(r_Aineq,1) = length(genes)-1;
        end
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*current_G(i) - 1)) = 1;
        end
        Aineq(r_Aineq,rNum + r) = length(current_G);
        bineq(r_Aineq,1) = 2*length(current_G)-1;
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*current_G(i) - 1)) = -1;
        end
        Aineq(r_Aineq,rNum + r) = -1 * length(current_G);
        bineq(r_Aineq,1) = -1 * length(current_G);
    end
    
    % up-regulation of complex isoenzymes
    if ChkList(r) == 4
        complexes = ParseCmplx(char(GeneRules{r}));
        current_G = [];
        for i = 1:size(complexes,1)
            genes = complexes{i,1};
            cmplx_inx = 0;
            for j = 1:variable_end
                if sum(variable(genes,j)) == length(genes)
                    if sum(variable(setdiff(1:GeneNum,genes),j)) == 0
                        cmplx_inx = j;
                        break
                    end
                end
            end
            if cmplx_inx == 0
                variable_end = variable_end + 1;
                variable(genes,variable_end) = 1;
                
                cmplx_inx = variable_end;
            end
            current_G = union(current_G,cmplx_inx);
            
            % 
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + GeneNum + genes(j)) = 1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*cmplx_inx) = -1 * length(genes);
            bineq(r_Aineq,1) = length(genes) - 1;
            
            % 
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + GeneNum + genes(j)) = -1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*cmplx_inx) = length(genes);
            bineq(r_Aineq,1) = 0;
        end
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*current_G(i)) = 1;
        end
        Aineq(r_Aineq,2*rNum + r) = length(current_G);
        bineq(r_Aineq,1) = length(current_G);
        
        % 
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*current_G(i)) = -1;
        end
        Aineq(r_Aineq,2*rNum + r) = -1 * length(current_G);
        bineq(r_Aineq,1) = -1;
    end
    
    genes = [];
end


%% unique the rows

[Aineq_unique,inx_unique,~] = unique(Aineq,'rows');
bineq_unique = bineq(inx_unique);

%% equality matrix

Aeq = zeros(mNum, size(Aineq_unique,2));
beq = zeros(mNum, 1);

% equality 1
Aeq(1:mNum,1:rNum) = modelIrrev.S;

r_Aeq = size(Aeq,1);
for r = 1:length(GeneRules)
    
    % rules with single genes: equalities 13 1nd 14
    if ChkList(r) == 1
        genes = erase(char(GeneRules{r}),"(");
        genes = erase(genes,")");
        
        g = find(strcmp(genes,GeneNames));
        
        % 
        r_Aeq = r_Aeq + 1;
        Aeq(r_Aeq,3*rNum + g) = 1;
        Aeq(r_Aeq,rNum + r) = 1;
        beq(r_Aeq,1) = 1;
        
        
        %
        r_Aeq = r_Aeq + 1;
        Aeq(r_Aeq,3*rNum + GeneNum + g) = 1;
        Aeq(r_Aeq,2*rNum + r) = 1;
        beq(r_Aeq,1) = 1;
        
        
    end
end

%% for additional variables
if 3*rNum+2*GeneNum < size(Aeq,2)
    LB(3*rNum+2*GeneNum+1:size(Aeq,2)) = zeros(size(Aeq,2)-(3*rNum+2*GeneNum),1);
    UB(3*rNum+2*GeneNum+1:size(Aeq,2)) = ones(size(Aeq,2)-(3*rNum+2*GeneNum),1);
end


%% MILP
f = zeros(size(Aeq,2),1);

% minimizing the number of genes
f(3*rNum+1:3*rNum+2*GeneNum) = 1;

intcon = rNum+1:size(Aeq,2);
options = optimoptions('intlinprog','IntegerPreprocess','none','HeuristicsMaxNodes', 100,'MaxTime',345600,'Heuristics','round','BranchRule','maxpscost','LPMaxIterations',3e7,'IntegerPreprocess', 'none','LPPreprocess','none');

[x,fval,flag] = intlinprog(f,intcon,Aineq_unique,bineq_unique,Aeq,beq,LB,UB,options);

DISP = [' The optimum value is : ',num2str(fval)];
disp(DISP)

Genes = x(3*rNum+1:3*rNum+2*GeneNum);
Genes(abs(Genes) < 10^(-6)) = 0;
Genes(abs(Genes - 1) < 10^(-6)) = 1;

Y = x(rNum+1:3*rNum);
Y(abs(Y) < 10^(-6)) = 0;
Y(abs(Y - 1) < 10^(-6)) = 1;

inx_modulate = mod(find(Y == 0),rNum);
inx_modulate(inx_modulate == 0) = rNum;
inx_modulate_type = floor(find(Y == 0)/rNum);

FinalRules = modelIrrev.rules(inx_modulate);


inx_modulate (all (cellfun ('isempty', FinalRules), 2), :) = [];  % remove empty cells
inx_modulate_type (all (cellfun ('isempty', FinalRules), 2), :) = [];  % remove empty cells
FinalRules (all (cellfun ('isempty', FinalRules), 2), :) = [];  % remove empty cells

% FinalRules gives the final list of reactions to be manipuated
FinalRules
%
Genes_inx = find(Genes > 0);
numel(Genes_inx);
numel(unique(mod(Genes_inx,GeneNum)));

FinalGeneList = GeneNames(unique(mod(Genes_inx,GeneNum)));
% FinalGeneList gives the final list of genes to be manipuated and
% gene_modulate_type shows the moduation type
FinalGeneList

FinalRules = ConvertEcoliRules(FinalRules,GeneNames);
gene_modulate_type = zeros(length(FinalGeneList),1);
for i = 1:length(FinalGeneList)
    for j = 1:length(FinalRules)
        if ~isempty(strfind(FinalRules{j},FinalGeneList{i}))
            gene_modulate_type(i) = inx_modulate_type(j);
            break
        end
    end
end

