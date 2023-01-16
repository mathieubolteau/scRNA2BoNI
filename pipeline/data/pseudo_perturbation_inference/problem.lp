%Generer combinaisons de n elements (n protéines)
{selprot(V):pert(E,V,S,C)} = k .

% Generer les vecteurs de expurbations correspondantes
aff(E,V,S,C) :- selprot(V), pert(E,V,S,C).

%Calculer les égalités entre affs de diffèrentes classe
egale(I,J,V) :- aff(I,V,S1,C1), aff(J,V,S2,C2), C1<C2, S1 = S2.

%Compter le nombre de egales entre 2 experiences de 2 classes differents, parmis les k protéines choisies
countegale(I,J,M) :- M={egale(I,J,_)}, aff(I,_,_,C1),aff(J,_,_,C2),C1<C2.

%Generer un prédicat affinite lorsque le nombre de coutnegale est egale à k pour des classes differentes.  Cela garantie que les affs soient egaux pour tous le protéines.
% Affinité entre les 2 exp si il y a k égalité en regardant les egale. On veut garder seulement les exp qui sont identiques entre les 2 classes et pour toutes les prot
0{affinite(I,J)}1 :- countegale(I,J,k),I!=J,aff(I,_,_,C1),aff(J,_,_,C2),C1<C2.

%Calculer les differences entre affs d'une meme classe
% Repérer les exp pour lesquelles une protéine (parmis celle selectionnée au début) n'est pas exprimée de la même manière chez les 2 exp

diff(I,J,V) :- aff(I,V,S1,C1), aff(J,V,S2,C2), C1==C2, S1!=S2, I<J. % Old version
% diff(I1,I2,V) :- aff(I1,V,S1,C), aff(I2,V,S2,C), S1!=S2, I1<I2.     % My version (MB)

%Compter le nombre de difference entre une meme classe
countdiff(I1,I2,M) :- M={diff(I1,I2,_)},I1<I2,aff(I1,_,_,C1),aff(I2,_,_,C2),C1==C2.

%Interdire qu'il existe un predicat countdiff de taille 0 pour la meme classe. Il suffit de le faire par proteine
% On gère les doublons, cad qu'un même vecteur soit présent pour plusieurs exp dans une meme classe
:- countdiff(I1,I2,0), aff(I1,_,_,C1),aff(I2,_,_,C2),C1==C2, affinite(I1,_), affinite(I2,_), I1<I2. % voir I!=J ?
:- countdiff(I1,I2,0), aff(I1,_,_,C1),aff(I2,_,_,C2),C1==C2, affinite(_,I1), affinite(_,I2), I1<I2. % RM ?

%%%%%%%%%% V4 -- Version 2  %%%%%%%%%%
selinput(G) :- selprot(G), input(G).

% Count the number of input genes expressed at 1 for each affinity
nb_selinput_expr1_by_aff(C, N) :- N={pert(C,G,1,_) : selinput(G)}, affinite(C,_).
% Forbid the affinity where the number of input genes expressed at 1 is inferior to 1
:- affinite(C,_), nb_selinput_expr1_by_aff(C,N), N < 1.


#maximize{1, I : affinite(I,_)}.

nb_aff(N) :- N= #count{1,I,J : affinite(I,J)}.

#show selprot/1.
#show affinite/2.
% %#show countaffinite/2.
#show nb_aff/1.
