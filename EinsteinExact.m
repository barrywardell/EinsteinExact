(* ::Package:: *)

(* An exact solutions thorn *)

Get["KrancThorn`"];
Needs["Metrics`"];

SetEnhancedTimes[False];


(**************************************************************************************)
(* Tensors *)
(**************************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  beta, dtbeta, g, k, Jac, xx, XX, betap, dtbetap, G, K
}];

Map[AssertSymmetricDecreasing, 
{
  g[la,lb], k[la,lb], G[la,lb], K[la,lb]
}];


(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)


admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}},
   {"admbase::dtlapse", {dtalp}},
   {"admbase::dtshift", {dtbetax,dtbetay,dtbetaz}}};

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  Jac[ui,lj], xx[ui], XX[ui], X, Y, Z, betap[ui], dtbetap[ui], G[li,lj], K[li,lk]
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

k11=kxx; k21=kxy; k22=kyy; k31=kxz; k32=kyz; k33=kzz;
g11=gxx; g21=gxy; g22=gyy; g31=gxz; g32=gyz; g33=gzz;
beta1=betax; beta2=betay; beta3=betaz;
dtbeta1=dtbetax; dtbeta2=dtbetay; dtbeta3=dtbetaz;


(**************************************************************************************)
(* Rotation matrices *)
(**************************************************************************************)

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

(**************************************************************************************)
(* Some helper functions *)
(**************************************************************************************)
makeKrancFriendly[var_Symbol] := var;
makeKrancFriendly[var_[x___]] :=Symbol[ToString[var]<>( StringJoin@@(ToString/@{x}))];
makeKrancFriendly[Derivative[y___][var_][x___]] :=Symbol["d"<>( StringJoin@@(ToString/@{y}))<>ToString[var]<>( StringJoin@@(ToString/@{x}))];
makeKrancFriendly[x_] := x;

(**************************************************************************************)
(* Create a Kranc thorn for a named spacetime *)
(**************************************************************************************)
idThorn[spacetime_] :=
  Module[{m, fourMetric, invFourMetric, threeMetric, lapse, shift,
          extrinsicCurvature, calc, calculations, parameters, coords, spatialCoords,
          metricShorthands, shortVars, krancShortVars, dShorthands, extraShorthands,
          extendedKeywordParameters},
  Print["Generating thorn for ", spacetime];
  m = GetMetric[spacetime];
  coords = ("Coordinates" /. m) /. {x -> X, y -> Y, z -> Z};
  If[coords =!= {t, X, Y, Z},
    Throw["Error, only metrics in Cartesian coordinates are supported"];
  ];
  spatialCoords = coords[[2;;]];

  (* Deal with shorthands - create new equations for all shorthands and their derivatives *)
  metricShorthands = (("Shorthands" /. m) /. "Shorthands" -> {}) /. {x -> X, y -> Y, z -> Z};
  shortVars = Flatten[{metricShorthands[[All,1]], D[metricShorthands[[All,1]], #]& /@ coords}];
  krancShortVars = (# -> makeKrancFriendly[#]) & /@ shortVars;
  dShorthands = Simplify[Flatten[D[metricShorthands, #]& /@ coords]];
  extraShorthands = (Join[metricShorthands, dShorthands]//. krancShortVars)/. (0->0) -> Sequence[];

  (* Compute lapse, shift, three metric and extrinsic curvature *)
  fourMetric = ("Metric" /. m) /. {x -> X, y -> Y, z -> Z};
  invFourMetric = Simplify[Inverse[fourMetric]];
  lapse = Simplify[1/Sqrt[-invFourMetric[[1,1]]]];
  shift = Simplify[lapse^2 invFourMetric[[1, 2;;4]]];
  threeMetric = fourMetric[[2;;4, 2;;4]];
  extrinsicCurvature = Simplify[Table[- 1/(2 lapse) (D[threeMetric[[i,j]], t]
    - Sum[D[threeMetric[[i,j]], spatialCoords[[k]]] shift[[k]], {k, 3}]
    - Sum[threeMetric[[i, k]] D[shift[[k]], spatialCoords[[j]]], {k, 3}]
    - Sum[threeMetric[[k, j]] D[shift[[k]], spatialCoords[[i]]], {k, 3}]),
    {j, 3}, {i, 3}]];

  (* Replace any shorthands with Kranc-friendly versions *)
  lapse = lapse //. krancShortVars;
  shift = shift //. krancShortVars;
  threeMetric = threeMetric //. krancShortVars;
  extrinsicCurvature = extrinsicCurvature //. krancShortVars;

  (* Get any necessary spacetime parameters *)
  parameters = ("Parameters" /. m)/. "Parameters" -> {};

  extendedKeywordParameters =
    Table[{Name -> paramName, AllowedValues -> {spacetime}},
      {paramName, {"ADMBase::initial_data", "ADMBase::initial_lapse", "ADMBase::initial_shift",
                   "ADMBase::initial_dtlapse", "ADMBase::initial_dtshift",
                   "ADMBase::evolution_method"}}];

  calc[when_] := {
    Name -> spacetime <> "_" <> when,
    Switch[when,
      "initial", Schedule -> {"in ADMBase_InitialData"},
      "always",  Schedule -> {"at ANALYSIS"},
      _, Throw["Unrecognised scheduling keyword"]],

    ConditionalOnKeyword ->
      {Switch[when,
         "initial", "initial_data",
         "always", "evolution_method",
         _, Throw["Unrecognised scheduling keyword"]],
       spacetime},

    Shorthands -> Join[shorthands, extraShorthands[[All,1]]],
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

      (* Add any shorthand equations *)
      extraShorthands,

      (* Compute unrotated variables *)
      Table[G[i,j]-> threeMetric[[i, j]], {j, 3}, {i, j, 3}],
      Table[K[i,j] -> extrinsicCurvature[[i,j]], {j, 3}, {i, j, 3}],
      Table[betap[i] -> shift[[i]], {i, 3}],
      Table[dtbetap[i] -> D[shift[[i]], t], {i, 3}],

      (* Computed rotated versions *)
      g[li,lj] -> Jac[um,li] Jac[un,lj] G[lm,ln],
      k[li,lj] -> Jac[um,li] Jac[un,lj] K[lm,ln],

      (* TODO: set these only if initial_lapse etc are set to
         <spacetime>.  This is not supported by Kranc at the
         moment. *)
      alp -> lapse,
      beta[ui] -> Jac[lj,ui] betap[uj],
      dtalp -> D[alp, t],
      dtbeta[ui] -> Jac[lj,ui] dtbetap[uj]
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
    ExtendedKeywordParameters -> extendedKeywordParameters,
    InheritedImplementations -> {"admbase"}];
];

spacetimes = {"GaugeWave", "KerrSchild", "Minkowski", "ShiftedGaugeWave"};

idThorn /@ spacetimes;



