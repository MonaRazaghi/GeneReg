clear;
clc;

load iJR904
model = iJR904;

[modelIrrev,~,~,irrev2rev] = convertToIrreversible(model);

[mNum,rNum] = size(modelIrrev.S);

biomass = find(modelIrrev.c == 1); % 
Obj = 329; % Ethanol exchange reaction

per_biomass = 0.2;
f_Obj = 0.3;
CC = 0.01;
L = 4*rNum;

load v_min
load v_max
load v0_l
load v0_u

v_min = Mini;
v_max = Maxi;
v0_l = Mini_bio;
v0_u = Maxi_bio;

for i = 1:length(v_min)
    if v0_l(i,1) < v_min(i,1) 
        v0_l(i,1) = v_min(i,1);
    end
end

for i = 1:length(v0_l)
    if v0_u(i,1) < v0_l(i,1) 
        v0_u(i,1) = v0_l(i,1);
    end
end

for i = 1:length(v0_u)
    if v_max(i,1) < v0_u(i,1) 
        v_max(i,1) = v0_u(i,1);
    end
end

modelIrrev.c = zeros(rNum,1);
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


GeneRules = modelIrrev.grRules;
ChkList = Rules_type(GeneRules); 

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


%% inequality matrix

Aineq = sparse((5*rNum + 2*num_rev + 1 + GeneNum), 3 * rNum + 2 * GeneNum);
bineq = zeros((5*rNum + 2*num_rev + 1 + GeneNum), 1);

% inequality 2
Aineq(1:rNum,1:rNum) = eye(rNum);
Aineq(1:rNum,rNum+1:2*rNum) = diag(dd-v_max);
bineq(1:rNum,1) = dd;

% inequality 1
Aineq(rNum+1:2*rNum,1:rNum) = -eye(rNum);
Aineq(rNum+1:2*rNum,rNum+1:2*rNum) = diag(v0_l-v_min);
bineq(rNum+1:2*rNum,1) = -v_min;

% inequality 4
Aineq(2*rNum+1:3*rNum,1:rNum) = eye(rNum);
Aineq(2*rNum+1:3*rNum,2*rNum+1:3*rNum) = diag(v_max-v0_u);
bineq(2*rNum+1:3*rNum,1) = v_max;

% inequality 3
Aineq(3*rNum+1:4*rNum,1:rNum) = -eye(rNum);
Aineq(3*rNum+1:4*rNum,2*rNum+1:3*rNum) = diag(v_min-uu);
bineq(3*rNum+1:4*rNum,1) = -uu;

% inequality 7
Aineq(4*rNum+1:5*rNum,rNum+1:2*rNum) = -eye(rNum);
Aineq(4*rNum+1:5*rNum,2*rNum+1:3*rNum) = -eye(rNum);
bineq(4*rNum+1:5*rNum,1) = -(diag(eye(rNum)));

% inequality 8
Aineq(5*rNum+1,rNum+1:3*rNum) = -1;
bineq(5*rNum+1,1) = L - 2*rNum;

% inequality 9 and 10
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

% genes
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

LB(biomass) = per_biomass * v0_u(biomass); % biomass
LB(Obj) = f_Obj * v_max(Obj); %  Obj

UB(1:rNum,1) = modelIrrev.ub;
UB(rNum+1:3*rNum+2*GeneNum,1) = ones(2*rNum+2*GeneNum,1);

%% relating genes to reactions by linear constraints

r_Aineq = size(Aineq,1); % last row index
c_Aineq = size(Aineq,2); % last col index
inx_LB = size(LB,1);

max_complexes = 100;
variable = zeros(GeneNum,max_complexes);
variable_end = 0;

for r = 1:length(GeneRules)
    
    % rules with or
    if ChkList(r) == 2
        t = erase(char(GeneRules{r}),"(");
        t = erase(t,")");
        t = strtrim(t);
        
        genes = extractingGenes2(t);
        
        % 5
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = 1;
        end
        Aineq(r_Aineq,rNum + r) = length(genes);
        bineq(r_Aineq,1) = 2*length(genes)-1;
        
        % 6
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = -1;
        end
        Aineq(r_Aineq,rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1 * length(genes);
        
        % 7
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = 1;
        end
        Aineq(r_Aineq,2*rNum + r) = length(genes);
        bineq(r_Aineq,1) = length(genes);
        
        % 8
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = -1;
        end
        Aineq(r_Aineq,2*rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1;
        
    end
    
    % rules with and
    if ChkList(r) == 3
        t = erase(char(GeneRules{r}),"(");
        t = erase(t,")");
        t = strtrim(t);
        genes = extractingGenes2(t);
        
        % 9
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = 1;
        end
        Aineq(r_Aineq,rNum + r) = length(genes);
        bineq(r_Aineq,1) = length(genes);
        
        % 10
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + g) = -1;
        end
        Aineq(r_Aineq,rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1;
        
        % 11
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = 1;
        end
        Aineq(r_Aineq,2*rNum + r) = length(genes);
        bineq(r_Aineq,1) = 2*length(genes)-1;
        
        % 12
        r_Aineq = r_Aineq + 1;
        for i = 1:length(genes)
            g = find(strcmp(genes(i),GeneNames));
            Aineq(r_Aineq,3*rNum + GeneNum + g) = -1;
        end
        Aineq(r_Aineq,2*rNum + r) = -1 * length(genes);
        bineq(r_Aineq,1) = -1 * length(genes);
        
    end
    
    % rules with conflicts down
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
            
            % 13
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + genes(j)) = 1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*cmplx_inx - 1)) = -1 * length(genes);
            bineq(r_Aineq,1) = 0;
            
            % 14
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + genes(j)) = -1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*cmplx_inx - 1)) = length(genes);
            bineq(r_Aineq,1) = length(genes)-1;
        end
        
        % 15
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*current_G(i) - 1)) = 1;
        end
        Aineq(r_Aineq,rNum + r) = length(current_G);
        bineq(r_Aineq,1) = 2*length(current_G)-1;
        
        % 16
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + (2*current_G(i) - 1)) = -1;
        end
        Aineq(r_Aineq,rNum + r) = -1 * length(current_G);
        bineq(r_Aineq,1) = -1 * length(current_G);
    end
    
    % rules with conflicts up
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
            
            % 17
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + GeneNum + genes(j)) = 1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*cmplx_inx) = -1 * length(genes);
            bineq(r_Aineq,1) = length(genes) - 1;
            
            % 18
            r_Aineq = r_Aineq + 1;
            for j = 1:length(genes)
                Aineq(r_Aineq,3*rNum + GeneNum + genes(j)) = -1;
            end
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*cmplx_inx) = length(genes);
            bineq(r_Aineq,1) = 0;
        end
        
        % 19
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*current_G(i)) = 1;
        end
        Aineq(r_Aineq,2*rNum + r) = length(current_G);
        bineq(r_Aineq,1) = length(current_G);
        
        % 20
        r_Aineq = r_Aineq + 1;
        for i = 1:length(current_G)
            Aineq(r_Aineq,3*rNum + 2*GeneNum + 2*current_G(i)) = -1;
        end
        Aineq(r_Aineq,2*rNum + r) = -1 * length(current_G);
        bineq(r_Aineq,1) = -1;
    end
    
    genes = [];
end

size(Aineq,1)

%% unique the rows

[Aineq_unique,inx_unique,~] = unique(Aineq,'rows');
bineq_unique = bineq(inx_unique);

size(Aineq_unique,1)
%% equality matrix

Aeq = zeros(mNum, size(Aineq_unique,2));
beq = zeros(mNum, 1);

Aeq(1:mNum,1:rNum) = modelIrrev.S;

r_Aeq = size(Aeq,1);
for r = 1:length(GeneRules)
    
    % rules with single genes
    if ChkList(r) == 1
        genes = erase(char(GeneRules{r}),"(");
        genes = erase(genes,")");
        
        g = find(strcmp(genes,GeneNames));
        
        % 1
        r_Aeq = r_Aeq + 1;
        Aeq(r_Aeq,3*rNum + g) = 1;
        Aeq(r_Aeq,rNum + r) = 1;
        beq(r_Aeq,1) = 1;
        
        %3
        r_Aeq = r_Aeq + 1;
        Aeq(r_Aeq,3*rNum + GeneNum + g) = 1;
        Aeq(r_Aeq,2*rNum + r) = 1;
        beq(r_Aeq,1) = 1;
        
    end
end

%%
if 3*rNum+2*GeneNum < size(Aeq,2)
    LB(3*rNum+2*GeneNum+1:size(Aeq,2)) = zeros(size(Aeq,2)-(3*rNum+2*GeneNum),1);
    UB(3*rNum+2*GeneNum+1:size(Aeq,2)) = ones(size(Aeq,2)-(3*rNum+2*GeneNum),1);
end


%% MILP
f = zeros(size(Aeq,2),1);
f(3*rNum+1:3*rNum+2*GeneNum) = 1;
intcon = rNum+1:size(Aeq,2);
options = optimoptions('intlinprog','IntegerPreprocess','none','HeuristicsMaxNodes', 100,'MaxTime',345600,'Heuristics','round','BranchRule','maxpscost','LPMaxIterations',3e5,'IntegerPreprocess', 'none','LPPreprocess','none');

[x,fval,~] = intlinprog(f,intcon,Aineq_unique,bineq_unique,Aeq,beq,LB,UB,options);


DISP = [' The optimum value is : ',num2str(fval)];
disp(DISP)

%% modifications
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

Genes_inx = find(Genes > 0);
numel(Genes_inx);
numel(unique(mod(Genes_inx,GeneNum)));

FinalGeneList = GeneNames(unique(mod(Genes_inx,GeneNum)));

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

    