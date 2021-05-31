g = DiagonalMatrix[{-1, a[t]^2/(1 - k*r^2), a[t]^2*r^2, 
    a[t]^2*r^2*Sin[\[Theta]]^2}];
xx = {t, r, \[Theta], \[Phi]};
U = {1, 0, 0, 0};
U1 = Table[
   Sum[g[[\[Mu], \[Nu]]]*U[[\[Nu]]], {\[Nu], 1, 4}], {\[Mu], 1, 4}];
T = Table[
   Table[\[Rho][t]*U1[[\[Mu]]]*U1[[\[Nu]]] + 
     p[t]*(g[[\[Mu], \[Nu]]] + U1[[\[Mu]]]*U1[[\[Nu]]]), {\[Nu], 1, 
     4}], {\[Mu], 1, 4}];
DifGeo[g, xx];

(*----------------------------------------------------------*)

g // MatrixForm
Ricci // MatrixForm
T // MatrixForm
EinsteinG // MatrixForm

(*----------------------------------------------------------*)

FullSimplify@
  Table[8 \[Pi]*Diagonal[T][[i]] == Diagonal[EinsteinG][[i]], {i, 1, 
    4}] // MatrixForm
