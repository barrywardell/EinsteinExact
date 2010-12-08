
Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];


(**************************************************************************************)
(* Tensors *)
(**************************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  h, hInv, phi, A, K, alpha, Gam, beta, betat, R, Rphi, gamma, bssnmom,  bssnCcons,
  g, k, AInv, DDphi, B, DDalpha, divA, gInv, divh,
  gPhys, kPhys, n, dir, dampG, dampDD, dampDDD, dampDG, ell,
  Lambda00, Lambda0, Lambda, g00, g0,
  Jac, xx, XX, XX2, X, G, K
}];

Map[AssertSymmetricDecreasing, 
{
  h[la,lb], hInv[ua,ub], A[la,lb], R[la,lb], Rphi[la,lb],
  g[la,lb], k[la,lb]
}];


(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)


admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

declaredGroups = Join[];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups];

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  ell[li], ss, rr, Lambda00, Lambda0[li], Lambda[ui,lj], g0[li], g00
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  {Name -> amp, Default -> 10^-1},
  {Name -> period, Default -> 1},
  {Name -> theta, Default -> 0},
  {Name -> phi, Default -> 0},
  {Name -> psi, Default -> 0},
  vx
};

intParameters = 
{

}

k11=kxx; k21=kxy; k22=kyy; k31=kxz; k32=kyz; k33=kzz;
g11=gxx; g21=gxy; g22=gyy; g31=gxz; g32=gyz; g33=gzz;
beta1=betax; beta2=betay; beta3=betaz;

Rphi =
 {{Cos[phi], Sin[phi], 0},
  {-Sin[phi], Cos[phi], 0},
  {0, 0, 1}};

Rtheta = 
  {{1, 0, 0}, 
   {0, Cos[theta], Sin[theta]}, 
   {0, -Sin[theta], Cos[theta]}};

Rpsi =
  {{Cos[psi], Sin[psi], 0},
  {-Sin[psi], Cos[psi], 0},
  {0, 0, 1}};

Rot = FullSimplify[Rpsi.Rtheta.Rphi];

gaugeWaveCalc[when_] := 
{
  Name -> "gauge_wave_" <> when,
  Switch[when,
    "initial", Schedule -> {"in ADMBase_PostInitial"},
    "always",  Schedule -> {"at ANALYSIS"},
    _, Throw["Unrecognised scheduling keyword"]],

  ConditionalOnKeyword -> {"when", when},  
  Shorthands -> {Jac[ui,lj], xx[ui], XX[ui], XX2, X, G[li,lj], K[li,lk]},
  Equations -> 
  {
    xx[1] -> x,
    xx[2] -> y,
    xx[3] -> z,

    Sequence@@Flatten[Table[Jac[i,j] -> Rot[[i,j]], {i,1,3}, {j,1,3}]],

    XX[ui] -> Jac[lj,ui] xx[uj],
    X -> XX[1],

    G[li,lj] -> Euc[li,lj],
    G11 -> (1-amp Sin[2 Pi (X - t)/period]),
    K[li,lj] -> 0,
    K11 -> -Pi amp/period Cos[2 Pi (X - t)/period]/Sqrt[1- amp Sin[2Pi(X-t)/period]],
    alp -> Sqrt[(1-amp Sin[2 Pi (X - t)/period])],
    beta[ui] -> 0,

    g[li,lj] -> Jac[li,um] Jac[lj,un] G[lm,ln],
    k[li,lj] -> Jac[li,um] Jac[lj,un] K[lm,ln]
  }
};

Print[gaugeWaveCalc["initial"]];

calculations = 
{
  gaugeWaveCalc["initial"],
  gaugeWaveCalc["always"]
};

keywordParameters = 
{
  {Name -> "when",
   AllowedValues -> {"initial", "always"},
   Default -> "initial"}
};

CreateKrancThornTT[groups, ".", "GaugeWave", 
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  IntParameters -> intParameters,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters,
  InheritedImplementations -> {"admbase"}];
