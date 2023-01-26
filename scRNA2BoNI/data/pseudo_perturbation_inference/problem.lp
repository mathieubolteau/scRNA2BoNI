% Search of the maximal set of pairs of cells where the expression 
% of the k genes selected are the same for the two cells of the pair.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I,J): a pair of cells where their expression for the k selected genes 
%     are the same for both cells. I and J belong to two different classes.


% Generate combinations of k genes
{selgene(G) : pert(C,G,S,CL)} = k .

% Generate the corresponding perturbation vectors
selpert(E,V,S,C) :- selgene(V), pert(E,V,S,C).

% Generate an input predicate for genes that are inputs
selinput(G) :- selgene(G), input(G).

% Generate a equal/3 predicate for cells I and J, from two different classes, where their expression for the gene G is equal
equal(I,J,G) :- selpert(I,G,S1,C1), selpert(J,G,S2,C2), C1<C2, S1==S2.

% Count the number of equal/3 predicates between 2 cells of different classes, among the k selected genes
countequal(I,J,M) :- M={equal(I,J,_)}, selpert(I,_,_,C1), selpert(J,_,_,C2), C1<C2.

% Generate, or not, an affinity/2 predicate when the number of countequal/3
% is equal to k for cells I and J of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I,J)}1 :- countequal(I,J,k), I!=J, selpert(I,_,_,C1), selpert(J,_,_,C2), C1<C2.

% Count the number of input genes expressed at 1 for each affinity/2
nbInputOnes(C,N) :- N={pert(C,G,1,_) : selinput(G)}, affinity(C,_).

% Forbid the affinity/2 where the number of input genes expressed at 1 is less than 1
:- affinity(C,_), nbInputOnes(C,N), N < 1.

% Generate a diff/3 predicate for cells I1 and I2, from a same class, where their expression for the gene G is different
diff(I1,I2,G) :- selpert(I1,G,S1,C1), selpert(I2,G,S2,C2), C1==C2, S1!=S2, I1<I2.

% Count the number of diff/3 predicates between 2 cells of the same class, among the k selected genes
countdiff(I1,I2,M) :- M={diff(I1,I2,_)}, I1<I2, selpert(I1,_,_,C1), selpert(I2,_,_,C2),C1==C2.

% Forbids countdiff/3 predicates for cells I1 and I2, from the first class
% where there is no difference in expression values for the k selected genes
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(I1,_), affinity(I2,_), I1<I2.

% Make the same for the second class
:- countdiff(I1,I2,0), selpert(I1,_,_,C1), selpert(I2,_,_,C2), C1==C2, affinity(_,I1), affinity(_,I2), I1<I2.

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I,_)}.

% Show selgene/1 and affinity/2 predicates
#show selgene/1.
#show affinity/2.
