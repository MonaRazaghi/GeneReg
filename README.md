# Feasibility of predicted metabolic engineering strategies

GeneReg is constraint-based approach that facilitates the design of feasible metabolic engineering strategies at the gene level and is readily applicable to large-scale metabolic networks.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

To run this code on your local machine, you need to have Matlab2017b. Also, the Optimization Toolbox and the COBRA toolbox v3.0 are needed to be installed.
The main code to run is called GeneReg.m, which is located in https://github.com/MonaRazaghi/GeneReg. Other .m files are necessary to be downloaded in the same path. 

## Running the tests

ecoli_core_model.mat is a small e.coli model provided in https://github.com/MonaRazaghi/GeneReg to run the code.
Some parameters that are mentioned in the manuscript can be set(changed) in the lines 14-20.
Here is the explanation on these parameters:

"biomass" refers to the index of biomass reaction in the model. In our sample this index is equal to 13.
"Obj" refers to the index of objective reaction in the model. In our sample to increase the production of succinate this index is set to 39.
"per_biomass" shows the factor of optimal  growth we want to have in our solution. Here it is set to 0.2.
"f_Obj" shows the factor that the flux towards the objective is at least this given factor of the maximum product yield. The default value of this parameter is set to 0.3.
"CC" represents the tunable parameter and it is set to 0.01 here.
"L" shows the limit on the total number of reaction manipulations. Here by setting this parameter to 4*(number of reactions), there will be no constraint on this parameter.


After running the GeneReg.m file, the code asks for the input model. The input model has to be compatible by COBRA toolbox. The sample model of ecoli_core_model.mat is compatible with the COBRA toolbox.
This model has 72 metabolites and 95 reactions, and for a model of this size it takes approximately five minutes to finish the optimization and to find the solution. In the output the optimal number of genes
need to be manipulated, the GPR rules for final reactions to manipulate and the final list of genes to manipulate are preseneted.
Moreover, the following variables provide some information on the results:

"FinalGeneList": the final gene list to manipulate
"gene_modulate_type": a vector of size FinalGeneList, which shows the manipulation type for each gene (0: downregulation, 1: upregulation)

"FinalRules": the GPR rule for final reactions to manipulate
"inx_modulate_type": a vector of size FinalRules, which shows the manipulation type for each reaction (0: downregulation, 1: upregulation)
"inx_modulate": a vector of size FinalRules, which shows the index of reactions in the model to manipulate (0: downregulation, 1: upregulation)

If there is no feasible solution for the optimization all above variables will be empty. 

## Authors: Zahra Razaghi-Moghadam, Zoran Nikoloski

