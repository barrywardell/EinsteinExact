(* ::Package:: *)

(* A gauge wave exact solution thorn *)

Get["KrancThorn`"];
Needs["Metrics`"];

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
  Jac, xx, XX, X, Y, Z, G, K
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
  {Name -> theta, Default -> 0},
  {Name -> phi, Default -> 0},
  {Name -> psi, Default -> 0},
  vx
};

intParameters =
{

};

keywordParameters =
{
  {Name -> "when",
   AllowedValues -> {"initial", "always", "dummy"},
   Default -> "initial"}
};

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

idThorn[spacetime_] :=
  Module[{m, fourMetric, threeMetric, lapse, shift, extrinsicCurvature, calc, calculations, parameters},
  m = GetMetric[spacetime];
  If[("Coordinates" /. m) =!= {t, x, y, z},
    Throw["Error, only metrics in Cartesian coordinates are supported"];
  ];

  parameters = ("Parameters" /. m)/. "Parameters" -> {};
  fourMetric = ("Metric" /. m) /. {x -> X, y -> Y, z -> Z};

  threeMetric = fourMetric[[2;;4, 2;;4]];
  lapse = -fourMetric[[1,1]];
  shift = fourMetric[[1, 2;;4]];
  extrinsicCurvature = - 1/2 D[threeMetric, t] / Sqrt[lapse];

  calc[when_] := {
    Name -> "gauge_wave_" <> when,
    Switch[when,
      "initial", Schedule -> {"in ADMBase_PostInitial"},
      "always",  Schedule -> {"at ANALYSIS"},
      _, Throw["Unrecognised scheduling keyword"]],

    ConditionalOnKeyword -> {"when", when},
    Shorthands -> {Jac[ui,lj], xx[ui], XX[ui], X, Y, Z, G[li,lj], K[li,lk]},
    Equations -> Flatten@
    {
      xx[1] -> x,
      xx[2] -> y,
      xx[3] -> z,

      Table[Jac[i,j] -> Rot[[i,j]], {i,1,3}, {j,1,3}],

      XX[ui] -> Jac[ui,lj] xx[uj],
      X -> XX[1],
      Y -> XX[2],
      Z -> XX[3],

      Table[Symbol["G"<>ToString[i]<>ToString[j]] -> threeMetric[[i, j]], {j, 3}, {i, 3}],
      Table[Symbol["K"<>ToString[i]<>ToString[j]] -> extrinsicCurvature[[i, j]], {j, 3}, {i, 3}],
      alp -> lapse,
      Table[Symbol["beta"<>ToString[i]] -> shift[[i]], {i, 3}],

      g[li,lj] -> Jac[um,li] Jac[un,lj] G[lm,ln],
      k[li,lj] -> Jac[um,li] Jac[un,lj] K[lm,ln]
    }
  };

  calculations =
  {
    calc["initial"],
    calc["always"]
  };

  CreateKrancThornTT[groups, "thorns", "InitialData_"<>spacetime,
    Calculations -> calculations,
    DeclaredGroups -> declaredGroupNames,
    IntParameters -> intParameters,
    RealParameters -> Join[realParameters, parameters],
    KeywordParameters -> keywordParameters,
    InheritedImplementations -> {"admbase"}];
];

spacetimes = {"GaugeWave", "KerrSchild", "Minkowski"};

idThorn /@ spacetimes;
