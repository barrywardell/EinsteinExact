(* ::Package:: *)

(* A gauge wave exact solution thorn *)

Get["KrancThorn`"];
Needs["Metrics`"];

SetEnhancedTimes[False];


(**************************************************************************************)
(* Tensors *)
(**************************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  beta, betat, g, k, Jac, xx, XX, G, K
}];

Map[AssertSymmetricDecreasing, 
{
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

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  Jac[ui,lj], xx[ui], XX[ui], X, Y, Z, G[li,lj], K[li,lk]
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters =
{
  {Name -> theta, Default -> 0},
  {Name -> phi, Default -> 0},
  {Name -> psi, Default -> 0}
};

keywordParameters =
{
  {Name -> "when",
   AllowedValues -> {"initial", "always"},
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
  lapse = -fourMetric[[1,1]]; (* from inverse metric *)
  shift = fourMetric[[1, 2;;4]]; (* from inverse metric *)
  extrinsicCurvature = - 1/2 D[threeMetric, t] / Sqrt[lapse]; (* shift terms *)

  calc[when_] := {
    Name -> spacetime <> when, (* rename *)
    Switch[when,
      "initial", Schedule -> {"in ADMBase_PostInitial"},
      "always",  Schedule -> {"at ANALYSIS"},
      _, Throw["Unrecognised scheduling keyword"]],

    ConditionalOnKeyword -> {"when", when},
    Shorthands -> shorthands,
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
      k[li,lj] -> Jac[um,li] Jac[un,lj] K[lm,ln],
      beta[ui] -> Jac[lj,ui] betap[uj]
      (* dtlapse and dtshift *)
    }
  };

  calculations =
  {
    calc["initial"],
    calc["always"]
  };

  CreateKrancThornTT[admGroups, "thorns", spacetime,
    Calculations -> calculations,
    RealParameters -> Join[realParameters, parameters],
    KeywordParameters -> keywordParameters,
    InheritedImplementations -> {"admbase"}];
];

spacetimes = {"GaugeWave", "KerrSchild", "Minkowski"};

idThorn /@ spacetimes;
