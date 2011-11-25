(* ::Package:: *)

(* EinsteinExact: An exact solutions thorn 
   Copyright 2011 Barry Wardell and Ian Hinder
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*)

Get["KrancThorn`"];

PrependTo[$Path, "."];
Needs["Metrics`"];

SetEnhancedTimes[False];


(**************************************************************************************)
(* Tensors *)
(**************************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  beta, dtbeta, g, k, Jac, InvJac, xx, XX, betap, dtbetap, G, K
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
  Jac[ui,lj], InvJac[ui,lj], xx[ui], XX[ui], X, Y, Z, betap[ui], dtbetap[ui], G[li,lj], K[li,lk]
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
  Module[{coordRule, coords, spatialCoords, fourMetric, invFourMetric,
          shorthandEquations, shorthandVars, dShorthands, simplifyhints, tf, cf, simpopts,
          threeMetric, lapse, shift, extrinsicCurvature, dtlapse, dtshift,
          krancShortVars, kranctf, parameters, extendedKeywordParameters, calc,
          calculations},

  Print["Generating thorn for ", spacetime];

  (* Load the spacetime: coordinates, metric, inverse metric *)
  coordRule = {x -> X, y -> Y, z -> Z, None -> {}};

  coords = MetricProperty[spacetime, "Coordinates"] /. coordRule;
  If[coords =!= {t, X, Y, Z},
    Throw["Error, only metrics in Cartesian coordinates are supported"];
  ];
  spatialCoords = coords[[2;;]];

  fourMetric = MetricProperty[spacetime, "Metric"]/. coordRule;
  invFourMetric = MetricProperty[spacetime, "InverseMetric"] /. coordRule;

  (* Create new equations for all shorthands and get rules for their derivatives *)
  shorthandEquations = MetricProperty[spacetime, "Shorthands"] /. coordRule;
  shorthandVars = shorthandEquations[[All,1]];

  (* Use simplification hints if provided by the metric to create TransformationFunctions for Simplify *)
  simplifyhints = MetricProperty[spacetime, "SimplifyHints"] /. coordRule;
  tf = {Automatic, Sequence@@(Function[{expr}, expr /. #] & /@
    Flatten[({simplifyhints, shorthandEquations, Reverse /@ shorthandEquations})])};

  (* Don't overcount the complexity of shorthand variables. They should really only be 1 leaf. *)
  cf = Function[{e}, LeafCount[e] + Apply[Plus, Map[Count[e, Blank[Head[#]], {0, Infinity}] (1 - LeafCount[#])&, shorthandVars]]];

  (* We pass both transformation and complexity functions to Simplify *)
  simpopts = Sequence[TransformationFunctions -> tf, ComplexityFunction -> cf];

  (* Get simplified expressions for the derivatives of the shorthands *)
  dShorthands = Simplify[Flatten[D[shorthandEquations, #] & /@ coords], simpopts] /. (0->0) -> Sequence[];

  (* If the inverse four metric is not provided then compute it *)
  If[invFourMetric === {},
    invFourMetric = Simplify[Inverse[fourMetric], simpopts];
  ];

  (* Compute lapse, shift, three metric and extrinsic curvature *)
  lapse = 1/Sqrt[-invFourMetric[[1, 1]]];
  shift = Simplify[lapse^2 invFourMetric[[1, 2 ;; 4]], simpopts];
  threeMetric = fourMetric[[2 ;; 4, 2 ;; 4]];
  extrinsicCurvature = -1/(2 alp) Simplify[Table[
     (D[threeMetric[[i, j]], t] -
       Sum[D[threeMetric[[i, j]], spatialCoords[[k]]] shift[[k]], {k, 3}] -
       Sum[threeMetric[[i, k]] D[shift[[k]], spatialCoords[[j]]], {k, 3}] -
       Sum[threeMetric[[k, j]] D[shift[[k]], spatialCoords[[i]]], {k, 3}]) /. dShorthands,
    {j, 3}, {i, 3}], simpopts];
  dtlapse = Simplify[D[lapse, t], simpopts];
  dtshift = Simplify[D[shift, t], simpopts];

  (* Create Kranc-friencly names for shorthands and corresponding transformation functions *)
  krancShortVars = (# -> makeKrancFriendly[#]) & /@ shorthandVars;
  kranctf = tf /. krancShortVars;

  (* Replace any shorthands with Kranc-friendly versions *)
  lapse = lapse /. krancShortVars;
  shift = shift /. krancShortVars;
  threeMetric = threeMetric /. krancShortVars;
  extrinsicCurvature = extrinsicCurvature /. krancShortVars;
  shorthandEquations = shorthandEquations /. krancShortVars;
  shorthandVars = shorthandVars /. krancShortVars;

  (* Get any necessary spacetime parameters *)
  parameters = MetricProperty[spacetime, "Parameters"] /. None -> {};

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

    Shorthands -> Join[shorthands, shorthandVars],
    Equations -> Flatten@
    {
      xx[1] -> x,
      xx[2] -> y,
      xx[3] -> z,

      (* The inverse of the rotation matrix is just its transpose *)
      Table[Jac[i,j] -> Rot[[i,j]], {i,1,3}, {j,1,3}],
      Table[InvJac[i,j] -> Jac[j,i], {i,1,3}, {j,1,3}],

      (* We compute the everything in the rotated coordinates *)
      XX[ui] -> Jac[ui,lj] xx[uj],
      X -> XX[1],
      Y -> XX[2],
      Z -> XX[3],

      (* Add any shorthand equations *)
      shorthandEquations,

      (* Compute rotated components *)
      (* TODO: set these only if initial_lapse etc are set to
         <spacetime>.  This is not supported by Kranc at the
         moment. *)
      alp -> lapse,
      dtalp -> dtlapse,
      Table[G[i,j]-> threeMetric[[i, j]], {j, 3}, {i, j, 3}],
      Table[K[i,j] -> extrinsicCurvature[[i,j]], {j, 3}, {i, j, 3}],
      Table[betap[i] -> shift[[i]], {i, 3}],
      Table[dtbetap[i] -> dtshift[[i]], {i, 3}],

      (* Rotate back to Cactus components *)
      g[li,lj] -> Jac[um,li] Jac[un,lj] G[lm,ln],
      k[li,lj] -> Jac[um,li] Jac[un,lj] K[lm,ln],
      beta[ui] -> InvJac[ui,lj] betap[uj],
      dtbeta[ui] -> InvJac[ui,lj] dtbetap[uj]
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



