DifGeo[g__, x__] := Module[
   {g0 = g, x0 = x, l = Length[x], invg = Inverse[g]}, 
   Chris = Table[
     Table[Table[
       Simplify@
        Sum[1/2*invg[[\[Sigma], \[Rho]]]*(D[g0[[\[Rho], \[Mu]]], 
             x0[[\[Nu]]]] + D[g0[[\[Nu], \[Rho]]], x0[[\[Mu]]]] - 
            D[g0[[\[Mu], \[Nu]]], x0[[\[Rho]]]]), {\[Rho], 1, 
          l}], {\[Nu], 1, l}], {\[Mu], 1, l}], {\[Sigma], 1, l}];
   (*Christoffel Symbol Subsuperscript[\[CapitalGamma],  \[Mu]\[Nu], \
\[Sigma]]=\[CapitalGamma][\[Sigma],\[Mu],\[Nu]]*)
   (*---------------------------------------------------------------------------*)
   
   Riemann = 
    Table[Table[
      Table[Table[
        Simplify[
         D[Chris[[\[Rho], \[Mu], \[Sigma]]], x0[[\[Nu]]]] - 
          D[Chris[[\[Rho], \[Nu], \[Sigma]]], x0[[\[Mu]]]] + 
          Sum[Chris[[\[Lambda], \[Sigma], \[Mu]]]*
             Chris[[\[Rho], \[Nu], \[Lambda]]] - 
            Chris[[\[Lambda], \[Sigma], \[Nu]]]*
             Chris[[\[Rho], \[Mu], \[Lambda]]], {\[Lambda], 1, 
            l}]], {\[Rho], 1, l}], {\[Sigma], 1, l}], {\[Nu], 1, 
       l}], {\[Mu], 1, l}];
   (*Riemann Curvature Tensor Subsuperscript[
   R, \[Mu]\[Nu]\[Sigma],    \[Rho]]=
   R[\[Mu]_,\[Nu]_,\[Sigma]_,\[Rho]_]*)
   (*---------------------------------------------------------------------------*)
   
   Ricci = 
    Table[Table[
      Sum[Riemann[[\[Mu], \[Nu], \[Sigma], \[Nu]]], {\[Nu], 1, 
        l}], {\[Sigma], 1, l}], {\[Mu], 1, l}];
   (*Ricci Tensor Subscript[R, \[Mu]\[Sigma]]=
   Ricci[\[Mu],\[Sigma]]*)
   (*---------------------------------------------------------------------------*)
   
   Rscalar = 
    Sum[Sum[invg[[\[Mu], \[Sigma]]]*Ricci[[\[Mu], \[Sigma]]], {\[Mu], 
       1, l}], {\[Sigma], 1, l}];
   (*Scalar Curvature*)
   (*---------------------------------------------------------------------------*)
   
   Weyl = 
    If[l < 3, Nothing, 
     Table[Table[
       Table[Table[
         Sum[Riemann[[\[Mu], \[Nu], \[Sigma], \[Eta]]]*
            g0[[\[Rho], \[Eta]]], {\[Eta], 1, l}] - 
          1/(n - 2)*(g0[[\[Mu], \[Sigma]]]*Ricci[[\[Rho], \[Nu]]] - 
             g0[[\[Mu], \[Rho]]]*Ricci[[\[Sigma], \[Nu]]] - 
             g0[[\[Nu], \[Sigma]]]*Ricci[[\[Rho], \[Mu]]] + 
             g0[[\[Nu], \[Rho]]]*Ricci[[\[Sigma], \[Mu]]]) + 
          1/((n - 1)*(n - 2))*
           Rscalar*(g0[[\[Mu], \[Sigma]]]*g0[[\[Rho], \[Nu]]] - 
             g0[[\[Mu], \[Rho]]]*g0[[\[Sigma], \[Nu]]]), {\[Rho], 1, 
          l}], {\[Sigma], 1, l}], {\[Nu], 1, l}], {\[Mu], 1, l}]];
   (*Weyl Tensor Subscript[C, \[Mu]\[Nu]\[Sigma]\[Rho]]=
   Weyl[\[Mu],\[Nu],\[Sigma],\[Rho]]*)
   (*---------------------------------------------------------------------------*)
   
   EinsteinG = FullSimplify[Ricci - 1/2*Rscalar*g0];
   (*Einstein Tensor Subscript[G, \[Mu]\[Nu]]=
   EinsteinG[\[Mu],\[Nu]]*)
   (*---------------------------------------------------------------------------*)
   
   killing = Table[Subscript[\[Xi], i] @@ x0, {i, 1, l}];
   killingeqn = {};
   prekillingeqn = 
    Table[Table[
      Sum[D[g0[[\[Nu], \[Sigma]]]*killing[[\[Sigma]]], x0[[\[Mu]]]] + 
        D[g0[[\[Mu], \[Sigma]]]*killing[[\[Sigma]]], x0[[\[Nu]]]] - 
        2 Chris[[\[Sigma], \[Mu], \[Nu]]]*
         Sum[g0[[\[Sigma], \[Rho]]]*killing[[\[Rho]]], {\[Rho], 1, 
           l}], {\[Sigma], 1, l}], {\[Mu], 1, l}], {\[Nu], 1, l}];
   Do[
    Do[
     If[MemberQ[killingeqn, prekillingeqn[[i, j]]], Nothing, 
      AppendTo[killingeqn, prekillingeqn[[i, j]]]],
     {i, 1, Dimensions[prekillingeqn][[1]]}], {j, 1, 
     Dimensions[prekillingeqn][[1]]}];
   (*Killing Equations*)
   (*---------------------------------------------------------------------------*)
   
   tran = Table[x0[[i]] -> x0[[i]][p], {i, 1, l}];
   geodesiceqn = 
    Table[D[x0[[\[Mu]]] /. tran, {p, 2}] + 
      Sum[Sum[(Chris[[\[Mu], \[Nu], \[Sigma]]] /. tran)*
         D[x0[[\[Nu]]] /. tran, p]*
         D[x0[[\[Sigma]]] /. tran, p], {\[Nu], 1, l}], {\[Sigma], 1, 
        l}], {\[Mu], 1, l}];
   (*Geodesic Equations*)
   ];
