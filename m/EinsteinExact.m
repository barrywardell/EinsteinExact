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
  beta, dtbeta, g, k, position, shiftadd, Jac, InvJac, xx, XX, betap, dtbetap, G, K
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
  position[ui], shiftadd[ui], Jac[ui,lj], InvJac[ui,lj], xx[ui], XX[ui], T, X, Y, Z, alpp, dtalpp, betap[ui], dtbetap[ui], G[li,lj], K[li,lk]
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters =
{
  (* Position at which the origin of the spacetime should be located *)
  {Name -> positiont, Default -> 0},
  {Name -> positionx, Default -> 0},
  {Name -> positiony, Default -> 0},
  {Name -> positionz, Default -> 0},
  (* Factor that should multiply the lapse *)
  {Name -> lapsefactor, Default -> 1},
  (* Shift that should be added *)
  {Name -> shiftaddx, Default -> 0},
  {Name -> shiftaddy, Default -> 0},
  {Name -> shiftaddz, Default -> 0},
  (* Angles by which the spacetime should be rotated *)
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
idThorn[spacetime_, thorn_] :=
  Module[{coordRule, coords, spatialCoords, fourMetric, invFourMetric,
          shorthandEquations, shorthandVars, dShorthands, simplifyhints, tf, cf, simpopts,
          threeMetric, lapse, shift, extrinsicCurvature, dtlapse, dtshift,
          krancShortVars, kranctf, parameters, extendedKeywordParameters, calc,
          calculations},

  Print["Generating thorn ", thorn, " for ", spacetime, " spacetime."];

  (* Load the spacetime: coordinates, metric, inverse metric *)
  coordRule = {t -> T, x -> X, y -> Y, z -> Z, None -> {}};

  coords = MetricProperty[spacetime, "Coordinates"] /. coordRule;
  If[coords =!= {T, X, Y, Z},
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

  (* We pass the transformation functions to Simplify *)
  simpopts = TransformationFunctions -> tf;

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
  extrinsicCurvature = -1/(2 lapse) Simplify[Table[
     (D[threeMetric[[i, j]], T] -
       Sum[D[threeMetric[[i, j]], spatialCoords[[k]]] shift[[k]], {k, 3}] -
       Sum[threeMetric[[i, k]] D[shift[[k]], spatialCoords[[j]]], {k, 3}] -
       Sum[threeMetric[[k, j]] D[shift[[k]], spatialCoords[[i]]], {k, 3}]) /. dShorthands,
    {j, 3}, {i, 3}], simpopts];
  dtlapse = Simplify[D[lapse, T] /. dShorthands, simpopts];
  dtshift = Simplify[D[shift, T] /. dShorthands, simpopts];

  (* Create Kranc-friencly names for shorthands and corresponding transformation functions *)
  krancShortVars = (# -> makeKrancFriendly[#]) & /@ shorthandVars;
  kranctf = tf /. krancShortVars;

  (* Replace any shorthands with Kranc-friendly versions *)
  lapse = lapse /. krancShortVars;
  shift = shift /. krancShortVars;
  dtlapse = dtlapse /. krancShortVars;
  dtshift = dtshift /. krancShortVars;
  threeMetric = threeMetric /. krancShortVars;
  extrinsicCurvature = extrinsicCurvature /. krancShortVars;
  shorthandEquations = shorthandEquations /. krancShortVars;
  shorthandVars = shorthandVars /. krancShortVars;

  (* Get any necessary spacetime parameters *)
  parameters = MetricProperty[spacetime, "Parameters"] /. None -> {};

  extendedKeywordParameters =
    Table[{Name -> paramName, AllowedValues -> {thorn}},
      {paramName, {"ADMBase::initial_data", "ADMBase::initial_lapse", "ADMBase::initial_shift",
                   "ADMBase::initial_dtlapse", "ADMBase::initial_dtshift",
                   "ADMBase::evolution_method"}}];

  calc[when_] := {
    Name -> thorn <> "_" <> when,
    Switch[when,
      "initial", Schedule -> {"in ADMBase_InitialData"},
      "always",  Schedule -> {"at CCTK_PRESTEP"},
      _, Throw["Unrecognised scheduling keyword"]],

    ConditionalOnKeyword ->
      {Switch[when,
         "initial", "initial_data",
         "always", "evolution_method",
         _, Throw["Unrecognised scheduling keyword"]],
       thorn},

    Shorthands -> Join[shorthands, shorthandVars],
    Equations -> Flatten@
    {
      xx[1] -> x,
      xx[2] -> y,
      xx[3] -> z,

      (* Position of origin *)
      position[1] -> positionx,
      position[2] -> positiony,
      position[3] -> positionz,

      (* Additional shift vector *)
      shiftadd[1] -> shiftaddx,
      shiftadd[2] -> shiftaddy,
      shiftadd[3] -> shiftaddz,

      (* The inverse of the rotation matrix is just its transpose *)
      Table[Jac[i,j] -> Rot[[i,j]], {i,1,3}, {j,1,3}],
      Table[InvJac[i,j] -> Jac[j,i], {i,1,3}, {j,1,3}],

      (* We compute everything in the rotated frame *)
      T -> lapsefactor (t - positiont),
      XX[ui] -> Jac[ui,lj] (xx[uj] - position[uj] - shiftadd[uj] T),
      X -> XX[1],
      Y -> XX[2],
      Z -> XX[3],

      (* Add any shorthand equations *)
      shorthandEquations,

      (* Compute rotated components *)
      (* TODO: set these only if initial_lapse etc are set to
         <spacetime>.  This is not supported by Kranc at the
         moment. *)
      Table[G[i,j]-> threeMetric[[i, j]], {j, 3}, {i, j, 3}],
      Table[K[i,j] -> extrinsicCurvature[[i,j]], {j, 3}, {i, j, 3}],
      alpp -> lapse,
      dtalpp -> dtlapse,
      Table[betap[i] -> shift[[i]], {i, 3}],
      Table[dtbetap[i] -> dtshift[[i]], {i, 3}],

      (* Rotate back to Cactus components *)
      g[li,lj] -> Jac[um,li] Jac[un,lj] G[lm,ln],
      k[li,lj] -> Jac[um,li] Jac[un,lj] K[lm,ln],
      alp -> lapsefactor alpp,
      dtalp -> lapsefactor dtalpp,
      beta[ui] -> InvJac[ui,lj] betap[uj] + shiftadd[ui],
      dtbeta[ui] -> InvJac[ui,lj] dtbetap[uj]
    }
  };

  calculations =
  {
    calc["initial"],
    calc["always"]
  };

  parameterConditions = 
  {{Parameter["initial_data"] != spacetime ||
    (Parameter["initial_lapse"] == spacetime &&
     Parameter["initial_shift"] == spacetime &&
     Parameter["initial_dtlapse"] == spacetime &&
     Parameter["initial_dtshift"] == spacetime),
   "When ADMBase::initial_data = \"" <> spacetime <> "\", the parameters ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse and ADMBase::initial_dtshift must also all be set to \"" <> spacetime <> "\""}};

  CreateKrancThornTT[admGroups, "../", thorn,
    Calculations -> calculations,
    RealParameters -> Join[realParameters, parameters],
    ExtendedKeywordParameters -> extendedKeywordParameters,
    ParameterConditions -> parameterConditions,
    InheritedImplementations -> {"admbase"},
    CSE -> True,
    UseLoopControl -> True,
    UseVectors -> True];
];

(* Create documentation stub *)
docGenerate[spacetime_, thorn_] :=
 Module[{shorthands, replaceShorthands, paramtext, shorthandtext, metrictextdoctext},
  parameters = MetricProperty[spacetime, "Parameters"] /. None -> {};
  paramtext =
    StringJoin @@ Riffle["\\code{" <> ToString[TeXForm[#]] <> "}" & /@ parameters, ", "];
  paramtext = StringReplace[paramtext, "\\epsilon" -> "epsilon"]; (* Horrible hack! *)

  shorthands = MetricProperty[spacetime, "Shorthands"] /. None -> {};
  replaceShorthands =
    Thread[shorthands[[All, 1]] ->
      Map[(# /. s_[___] :> s) &, shorthands[[All, 1]]]];
  shorthandtext =
    StringJoin[("\\begin{equation}\n" <> ToString[TeXForm[Equal @@ #]] <>
      "\n\\end{equation}\n") & /@ shorthands];

  metrictext =
    ToString[TeXForm[MetricProperty[spacetime, "Metric"] /. replaceShorthands]];

  doctext = "
\\subsection{" <> spacetime <> "}
{\\bf Description:}  " <>
ToString[TeXForm[MetricProperty[spacetime, "Description"]]] <> " \\\\
{\\bf Thorn name:} " <> thorn <>
If[paramtext =!= "", " \\\\
{\\bf Parameters:} " <> paramtext, ""] <> " \\\\
{\\bf Metric:} \\\\
\\begin{equation}
g_{ab} =
" <> metrictext <> "
\\end{equation}
" <> If[shorthandtext =!= "", "where\n" <> shorthandtext, ""];

  PutAppend[OutputForm[doctext], "../doc/spacetimes.tex"];
]

spacetimes = {"GaugeWave", "KerrSchild", "Minkowski", "ShiftedGaugeWave", "Vaidya", "ModifiedSchwarzschildBL"};
thorns     = {"GaugeWave", "KerrSchild", "Minkowski", "ShiftedGaugeWave", "Vaidya2", "ModifiedSchwarzschildBL"};

MapThread[idThorn, {spacetimes, thorns}];

DeleteFile["../doc/spacetimes.tex"];
MapThread[docGenerate, {spacetimes, thorns}];
