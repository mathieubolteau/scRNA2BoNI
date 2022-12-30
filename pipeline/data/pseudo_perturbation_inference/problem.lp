%#const k=3. % toy
%#const k=10. % larger

{selinput(G):input(G), selprot(G)}=i.

%Generer combinaisons de n elements (n protéines)
{selprot(V):pert(E,V,S,C)} = k .

% Generer tous les sous ensembles d'experiences possibles 
affi(E) :- pert(E,_,_,_).

% Generer les vecteurs de expurbations correspondantes
aff(E,V,S,C) :- selprot(V), pert(E,V,S,C).

%Eliminer les affs qui ne contienent pas l'ensemble de protéines qui m'interesse
:- aff(E,V,S,C), not selprot(V). % REDONDANT ?


%Calculer les égalités entre affs de diffèrentes classe
egale(I,J,V) :- aff(I,V,S1,C1), aff(J,V,S2,C2), C1<C2, S1 = S2.

%Compter le nombre de egales entre 2 experiences de 2 classes differents, parmis les k protéines choisies
countegale(I,J,M) :- M={egale(I,J,_)}, I!=J, affi(I), affi(J), aff(I,_,_,C1),aff(J,_,_,C2),C1<C2.

%Generer un prédicat affinite lorsque le nombre de coutnegale est egale à k pour des classes differentes.  Cela garantie que les affs soient egaux pour tous le protéines.
% Affinité entre les 2 exp si il y a k égalité en regardant les egale. On veut garder seulement les exp qui sont identiques entre les 2 classes et pour toutes les prot
0{affinite(I,J)}1 :- countegale(I,J,k),I!=J,affi(I),affi(J),aff(I,_,_,C1),aff(J,_,_,C2),C1<C2.

%Calculer les differences entre affs d'une meme classe
% Repérer les exp pour lesquelles une protéine (parmis celle selectionnée au début) n'est pas exprimée de la même manière chez les 2 exp

diff(I,J,V) :- aff(I,V,S1,C1), aff(J,V,S2,C2), C1==C2, S1!=S2, I<J. % Old version
% diff(I1,I2,V) :- aff(I1,V,S1,C), aff(I2,V,S2,C), S1!=S2, I1<I2.     % My version (MB)

%Compter le nombre de difference entre une meme classe
countdiff(I1,I2,M) :- M={diff(I1,I2,_)},I1<I2,affi(I1),affi(I2),aff(I1,_,_,C1),aff(I2,_,_,C2),C1==C2.

%Interdire qu'il existe un predicat countdiff de taille 0 pour la meme classe. Il suffit de le faire par proteine
% On gère les doublons, cad qu'un même vecteur soit présent pour plusieurs exp dans une meme classe
:- countdiff(I1,I2,0), affi(I1),affi(I2),aff(I1,_,_,C1),aff(I2,_,_,C2),C1==C2, affinite(I1,_), affinite(I2,_), I1<I2. % voir I!=J ?
:- countdiff(I1,I2,0), affi(I1),affi(I2),aff(I1,_,_,C1),aff(I2,_,_,C2),C1==C2, affinite(_,I1), affinite(_,I2), I1<I2. % RM ?



%%%%%%%%%% V2 %%%%%%%%ù

% Forbid vectors composed only of '0' => OK : verif : clingo -c k=4 instances/toy/falsified_toy_instance.lp instances/toy/toy_input_genes.lp V2_gene_cells_selection.lp
sumvector(I,J,N) :- N=#sum{S : aff(I,V,S,_)} , affinite(I,J).
:- sumvector(_,_,N), N<1.    
% Explorer N=0
% :- sumvector(_,_,0).
% :- sumvector(I,J,N), N<1, affinite(I,J). 



%%%%%%%%%% V4 -- Version 1  %%%%%%%%%%
% % Forbid model where no input is present in the selected genes

% selinput(G) :- selprot(G), input(G).
% :- not selinput(_).


% % Define selinput expressed to 1
% selinput_expr_1(G) :- selinput(G), pert(C,G,1,_), affinite(C,_). % Pour 1 affinite au moins (pas toutes, si on veut le faire: modifier cette ligne)

% % Mandatory at least 1 input expressed at 1
% nb_selinput_expr_1(N) :- N={selinput_expr_1(_)}.
% :- nb_selinput_expr_1(N), N==0.


%%%%%%%%%% V4 -- Version 2  %%%%%%%%%%
selinput(G) :- selprot(G), input(G).
:- not selinput(_). % OK => /home/e21g017n/Nextcloud/work/gitlab_repos/pseudo-perturbations/run_new_code/tests
                    % test_underscore2_50_l5.out

                    % To discuss: really necessary ??

% Count the number of input genes expressed at 1 for each affinity
nb_selinput_expr1_by_aff(C, N) :- N={pert(C,G,1,_) : selinput(G)}, affinite(C,_).
% Forbid the affinity where the number of input genes expressed at 1 is inferior to 1
:- affinite(C,_), nb_selinput_expr1_by_aff(C,N), N < 1.



%%%% Redondancies %%%

redondancyJ(I,J2) :- countegale(I,J2,k),I!=J1, I!=J2,affi(I),affi(J2),affi(J1), J1!=J2,aff(I,_,_,C1),aff(J2,_,_,C2),C1<C2, affinite(I,J1).
redondancyI(I2,J) :- countegale(I2,J,k),I1!=J, I2!=J,affi(I1),affi(I2),affi(J), I1!=I2,aff(I2,_,_,C1),aff(J,_,_,C2),C1<C2, affinite(I1,J).

% redondancyInferred(I2,J2) :- redondancyJ(I,J2), redondancyI(I2, J) , not affinite(I2,J2), countegale(I2,J2,k).

redondancy(I,J) :- redondancyI(I,J).
redondancy(I,J) :- redondancyJ(I,J).
% redondancy(I,J) :- redondancyInferred(I,J).

#maximize{1@1,I,J : redondancy(I,J)}.


nb_redondancyI(N) :- N= #count{1,I,J : redondancyI(I,J)}.
nb_redondancyJ(N) :- N= #count{1,I,J : redondancyJ(I,J)}.
% nb_redondancyInferred(N) :- N= #count{1,I,J : redondancyInferred(I,J)}.

nb_redondancy(N) :- N= #count{1,I,J : redondancy(I,J)}.

%le nombre maximale de aff
% #maximize{1@2, I : affinite(I,_)}.

nb_aff(N) :- N= #count{1,I,J : affinite(I,J)}.

% #show compteraffs/2.
% #show egale/3.
#show selprot/1.
#show affinite/2.
% %#show countaffinite/2.

#show nb_aff/1.

% #show selinput/1.
% #show nb_sel_input/1.

% #show nb_selinput_expr1_by_aff/2.


% #show sumvector/3.
% nb(N) :- N={sumvector(_,_,_)}.
% #show nb/1.

% #show redondancyI/1.
#show redondancyI/2.
#show redondancyJ/2.
% #show redondancyI/3.
% #show redondancyJ/3.

% #show redondancyInferred/2.
#show nb_redondancyI/1.
#show nb_redondancyJ/1.
% #show nb_redondancyInferred/1.

#show redondancy/2.
#show nb_redondancy/1.
% 

% selprot("AES"). selprot("CD44"). selprot("ERCC4"). selprot("SOD2"). selprot("EGF").
% selprot(gene1). selprot(gene7). selprot(gene4).


% selprot("AES"). selprot("CD44"). selprot("ERCC4"). selprot("SOD2"). selprot("IGF2").
% selprot("AES"). selprot("CD44"). selprot("ERCC4"). selprot("SOD2"). selprot("FOXA1").











% redondancyJ(J1,J2, I1) :- countdiff(J1,J2,0), affinite(I1, J1), J1!=J2.%, affi(I1), affi(I2), affi(J), aff(I1,_,_,C),aff(I2,_,_,C), aff(J,_,_,C2), C!=C2.
% redondancyI(I1,I2,J) : countdiff(I1,I2,0), affinite(I1, J), I1!=I2.
% % redondancyI(I,J2) :- countdiff(J1,J2,0), affi(J1), affi(J2), affi(J), aff(J1,_,_,C),aff(J2,_,_,C), aff(I,_,_,C2), C!=C2, affinite(I, J1), J1!=J2.

% % redondancyInferred(I2,J2) :- redondancyI(I1,I2,J), redondancyJ(J1,J2,I) , not affinite(I2,J2), countdiff(I1,J1,0).

% redondancy(I2,J) :- redondancyI(I1,I2,J) , affi(I1), affi(I2), affi(J).
% redondancy(I,J2) :- redondancyJ(J1,J2, I), affi(I), affi(J1), affi(J2).


% redondancy(I2,J) :- redondancyI(I1,I2), affinite(I1, J).%, aff(I1,_,_,C1),aff(I2,_,_,C1), aff(J,_,_,C2), C1!=C2, I1!=I2.
% redondancy(I2,J) :- redondancyI(I2,I1), affinite(I1, J).%, aff(I1,_,_,C1),aff(I2,_,_,C1), aff(J,_,_,C2), C1!=C2, I1!=I2.

% redondancy(I,J2) :- redondancyJ(J1,J2), affinite(I, J1).%, aff(J1,_,_,C2),aff(J2,_,_,C2), aff(I,_,_,C1), C1!=C2, J1!=J2.
% redondancy(I,J2) :- redondancyJ(J2,J1), affinite(I, J1).%, aff(J1,_,_,C2),aff(J2,_,_,C2), aff(I,_,_,C1), C1!=C2, J1!=J2.



% redondancyI(I,J2) :- countegale(I,J2,k),I!=J1, I!=J2,affi(I),affi(J2),affi(J1), J1!=J2,aff(I,_,_,C1),aff(J2,_,_,C2),C1<C2, redondancyI(I,J1).
% redondancyJ(I2,J) :- countegale(I2,J,k),I1!=J, I2!=J,affi(I1),affi(I2),affi(J), I1!=I2,aff(I2,_,_,C1),aff(J,_,_,C2),C1<C2, affiredondancyJnite(I1,J).

% nb_redondancyI(N) :- N= #count{1,I,J,K : redondancyI(I,J,K)}.
% nb_redondancyJ(N) :- N= #count{1,I,J,K : redondancyJ(I,J,K)}.
