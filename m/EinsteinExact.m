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
  alp, dtalp, betal, beta, dtbetal, dtbeta, dbetal, dbeta,
  g, dtg, dg, gu, dtgu, dgu, k
}];

AssertSymmetricIncreasing[g[la,lb], la, lb];
AssertSymmetricIncreasing[dtg[la,lb], la, lb];
AssertSymmetricIncreasing[dg[la,lb,lc], la, lb];
AssertSymmetricIncreasing[gu[ua,ub], ua, ub];
AssertSymmetricIncreasing[dtgu[ua,ub], ua, ub];
AssertSymmetricIncreasing[dgu[ua,ub,lc], ua, ub];
AssertSymmetricIncreasing[k[la,lb], la, lb];


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

(* Transform pseudo-tensors to symbols, e.g. g4[1,2] -> g412 *)
arrayToSymbolRules =
  Module[{names},
         names = {invXform1L, invXform2L, invXformL, dInvXformL,
                  xx, txx, tg4, tdg4, g4, dg4};
         {name_[inds__Integer] :>
          Symbol[StringJoin[Map[ToString, {name, inds}]]]
          /; MemberQ[names, name]}];
(* Handle symmetries explicitly *)
symmetryRules =
  Flatten[{Table[tg4[i,j]    -> tg4[j,i]   , {i, 0, 3}, {j, 0, i-1}],
           Table[tdg4[i,j,k] -> tdg4[j,i,k], {i, 0, 3}, {j, 0, i-1}, {k, 0, 3}],
           Table[g4[i,j]     -> g4[j,i]    , {i, 0, 3}, {j, 0, i-1}],
           Table[dg4[i,j,k]  -> dg4[j,i,k] , {i, 0, 3}, {j, 0, i-1}, {k, 0, 3}]
          }] /. arrayToSymbolRules;

shorthands =
  Join[
    {
      T, X, Y, Z,
      betal[li], detg, gu[ui,uj], betasq, dtg[li,lj], dg[li,lj,lk], dtgu[ui,uj],
      dgu[ui,uj,lk], dtbetal[li], dbetal[li,lj], dbeta[ui,lj], dtbetasq
    },
    Flatten[{
      Table[invXform1L[i,j], {i, 0, 3}, {j, 0, 3}],
      Table[invXform2L[i,j], {i, 0, 3}, {j, 0, 3}],
      Table[invXformL[i,j], {i, 0, 3}, {j, 0, 3}],
      Table[dInvXformL[i,j,k], {i, 0, 3}, {j, 0, 3}, {k, 0, 3}],
      Table[xx[i], {i, 0, 3}],
      Table[txx[i], {i, 0, 3}],
      Table[tg4[i,j], {i, 0, 3}, {j, i, 3}],
      Table[tdg4[i,j,k], {i, 0, 3}, {j, i, 3}, {k, 0, 3}],
      Table[g4[i,j], {i, 0, 3}, {j, i, 3}],
      Table[dg4[i,j,k], {i, 0, 3}, {j, i, 3}, {k, 0, 3}]
    }]] /. arrayToSymbolRules;

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters =
{
  (* Time offset *)
  {Name -> timeoffset, Default -> 0},
  (* Position of origin *)
  {Name -> positionx, Default -> 0},
  {Name -> positiony, Default -> 0},
  {Name -> positionz, Default -> 0},
  (* Rotation by Euler angles *)
  {Name -> theta, Default -> 0},
  {Name -> phi, Default -> 0},
  {Name -> psi, Default -> 0},
  (* Lapse factor *)
  {Name -> lapsefactor, Default -> 1},
  (* Lorentz boost *)
  {Name -> boostx, Default -> 0},
  {Name -> boosty, Default -> 0},
  {Name -> boostz, Default -> 0},
  (* Shift modification *)
  {Name -> shiftaddx, Default -> 0},
  {Name -> shiftaddy, Default -> 0},
  {Name -> shiftaddz, Default -> 0}
};

k11=kxx; k12=kxy; k22=kyy; k13=kxz; k23=kyz; k33=kzz;
g11=gxx; g12=gxy; g22=gyy; g13=gxz; g23=gyz; g33=gzz;
beta1=betax; beta2=betay; beta3=betaz;
dtbeta1=dtbetax; dtbeta2=dtbetay; dtbeta3=dtbetaz;



(******************************************************************************)
(* Transformation matrices *)
(******************************************************************************)

Print["Generating transformations"];

(* Could be Identity, Simplify, FullSimplify *)
xformSimplify = FullSimplify;

(* Rotation about Euler angles {theta, phi, psi} *)
xformEuler = Module[{delta4, Rphi, Rtheta, Rpsi, Rot, xform},
                    delta4 = IdentityMatrix[4];
                    Rphi = {{Cos[phi], -Sin[phi], 0},
                            {Sin[phi], Cos[phi], 0},
                            {0, 0, 1}};
                    Rtheta = {{1, 0, 0}, 
                              {0, Cos[theta], -Sin[theta]}, 
                              {0, Sin[theta], Cos[theta]}};
                    Rpsi = {{Cos[psi], -Sin[psi], 0},
                            {Sin[psi], Cos[psi], 0},
                            {0, 0, 1}};
                    Rot = Rphi.Rtheta.Rpsi;
                    xform = Table[If[i==0 || j==0, delta4[[i+1,j+1]],
                                     Rot[[i,j]]],
                                  {i,0,3}, {j,0,3}];
                    xformSimplify[xform]];

(* Modify lapse *)
xformSlowdown = Module[{delta4, xform},
                       delta4 = IdentityMatrix[4];
                       xform = Table[If[i==0 && j==0, lapsefactor,
                                        delta4[[i+1,j+1]]],
                                     {i,0,3}, {j,0,3}];
                       xformSimplify[xform]];

(* Not implemented: deformation, cushion, shear *)

(* Lorentz Boost *)
xformBoost = Module[{eps, delta3, delta4,
                     gamma, boost, xform0, xform1, xform3, xform},
                    delta3 = IdentityMatrix[3];
                    delta4 = IdentityMatrix[4];
                    boost = {boostx, boosty, boostz};
                    gamma = 1 / Sqrt[1 - Total[boost^2]];
                    xform0 = gamma;
                    xform1 = boost gamma;
                    xform3 = (delta3 +
                              gamma^2 / (gamma+1) Outer[Times, boost, boost]);
                    xform = Table[If[i==0,
                                     If[j==0, xform0, xform1[[j]]],
                                     If[j==0, xform1[[i]], xform3[[i,j]]]],
                                  {i,0,3}, {j,0,3}];
                    xformSimplify[xform]];

(* Modify shift *)
xformVelocity = Module[{delta4, velocity, xform1, xform},
                       delta4 = IdentityMatrix[4];
                       velocity = {shiftaddx, shiftaddy, shiftaddz};
                       xform1 = velocity;
                       xform = Table[If[i>0 && j==0, xform1[[i]],
                                        delta4[[i+1,j+1]]],
                                     {i,0,3}, {j,0,3}];
                       xformSimplify[xform]];

(* The order in which these transformations are chained does matter --
   the order below has been chosen carefully to "make sense" *)
(* We keep the transformation in two pieces that we multiply at run
   time to speed up code generation *)
xform1 = xformVelocity . xformBoost . xformSlowdown;
xform2 = xformEuler;
xform = xform1 . xform2;

(* Inverse transformation *)
invXform1 = Inverse[xform1];
invXform2 = Inverse[xform2];
invXform = invXform2 . invXform1;

(* Note: indices are xform[ua,lb] and invXform[ua,lb] *)

coordOffset = {timeoffset, positionx, positiony, positionz};
originalCoords = {t,x,y,z};

(* Derivative of inverse transformation *)
dInvXform = Table[D[invXform[[i, j]], originalCoords[[k]]],
                  {i, 4}, {j, 4}, {k, 4}];
(* This is currently always zero -- abort if not *)
If[Count[Flatten[dInvXform], Except[0]] > 0, Throw["Error, dInvXform is non-zero"]]

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
          dFourMetric,
          shorthandEquations, shorthandVars, dShorthands, simplifyhints, tf, cf, simpopts,
          (* lapse, shift, threeMetric, dtlapse, dtshift, dtthreeMetric,
          dxshift, dxthreeMetric, extrinsicCurvature, *)
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

  fourMetric = MetricProperty[spacetime, "Metric"] /. coordRule;

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

  (* Simplify an expression *)
  simp[expr_] := Simplify[expr /. dShorthands, simpopts];

  dFourMetric = simp[Table[D[fourMetric[[i, j]], coords[[k]]],
                           {i, 4}, {j, 4}, {k, 4}]];

  (* Create Kranc-friencly names for shorthands and corresponding transformation functions *)
  krancShortVars = (# -> makeKrancFriendly[#]) & /@ shorthandVars;
  kranctf = tf /. krancShortVars;

  (* Replace any shorthands with Kranc-friendly versions *)
  fourMetric = fourMetric /. krancShortVars;
  dFourMetric = dFourMetric /. krancShortVars;
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
    Equations -> (Flatten@
    {
      (* TODO: these transformations are the same for every grid point;
         do not recalculate them *)
      Table[invXform1L[i,j] -> invXform1[[i+1, j+1]], {i, 0, 3}, {j, 0, 3}],
      Table[invXform2L[i,j] -> invXform2[[i+1, j+1]], {i, 0, 3}, {j, 0, 3}],
      Table[invXformL[i,j] -> Sum[invXform2L[i,k] invXform1L[k,j], {k, 0, 3}],
            {i, 0, 3}, {j, 0, 3}],
      (* This is always zero
      Table[dInvXformL[i,j,k] -> dInvXform[[i+1, j+1, k+1]],
            {i, 0, 3}, {j, 0, 3}, {k, 0, 3}],
      *)

      Table[xx[i] -> originalCoords[[i+1]] - coordOffset[[i+1]], {i, 0, 3}],
      Table[txx[i] -> Sum[invXformL[i,j] xx[j], {j, 0, 3}], {i, 0, 3}],

      T -> txx[0],
      X -> txx[1],
      Y -> txx[2],
      Z -> txx[3],

      (* Add any shorthand equations *)
      shorthandEquations,

      (* Compute components *)
      Table[tg4[i,j] -> fourMetric[[i+1, j+1]], {i, 0, 3}, {j, i, 3}],
      Table[tdg4[i,j,k] -> dFourMetric[[i+1, j+1, k+1]],
            {i, 0, 3}, {j, i, 3}, {k, 0, 3}],
      Table[g4[i,j] -> Sum[invXformL[k,i]
                           Sum[invXformL[l,j] tg4[k,l], {l, 0, 3}],
                           {k, 0, 3}],
            {i, 0, 3}, {j, i, 3}],
      (* dInvXform is always zero *)
      Table[dg4[i,j,k] -> ((* + Sum[dInvXformL[l,i,k]
                                 Sum[invXformL[m,j] tg4[l,m], {m, 0, 3}],
                                 {l, 0, 3}] *)
                           (* + Sum[invXformL[l,i]
                                 Sum[dInvXformL[m,j,k] tg4[l,m], {m, 0, 3}],
                                 {l, 0, 3}] *)
                           + Sum[invXformL[l,i]
                                 Sum[invXformL[m,j] 
                                     Sum[invXformL[n,k] tdg4[l,m,n], {n, 0, 3}],
                                     {m, 0, 3}],
                                 {l, 0, 3}]),
            {i, 0, 3}, {j, i, 3}, {k, 0, 3}],

      (* Compute lapse, shift, three metric and extrinsic curvature *)
      (* TODO: set these only if initial_lapse etc are set to
         <spacetime>. This is not supported by Kranc at the moment. *)
      Table[betal[i] -> g4[0,i], {i, 3}],
      Table[g[i,j] -> g4[i,j], {i, 3}, {j, i, 3}],
      detg -> Det[MatrixOfComponents[g[li,lj]]],
      gu[ui,uj] -> (1/detg Det[MatrixOfComponents[g[li,lj]]]
                    MatrixInverse[g[ui,uj]]),
      beta[ui] -> gu[ui,uj] betal[lj],
      betasq -> betal[li] beta[ui],
      alp -> Sqrt[betasq - g4[0,0]],

      Table[dtg[i,j] -> dg4[i,j,0], {i, 3}, {j, i, 3}],
      Table[dg[i,j,k] -> dg4[i,j,k], {i, 3}, {j, i, 3}, {k, 3}],
      dtgu[ui,uj] -> - gu[ui,uk] gu[uj,ul] dtg[lk,ll],
      dgu[ui,uj,lk] -> - gu[ui,ul] gu[uj,um] dg[ll,lm,lk],
      Table[dtbetal[i] -> dg4[0,i,0], {i, 3}],
      Table[dbetal[i,j] -> dg4[0,i,j], {i, 3}, {j, 3}],
      dtbeta[ui] -> dtgu[ui,uj] betal[lj] + gu[ui,uj] dtbetal[lj],
      dbeta[ui,lj] -> dgu[ui,uk,lj] betal[lk] + gu[ui,uk] dbetal[lk,lj],
      dtbetasq -> dtbetal[li] beta[ui] + betal[li] dtbeta[ui],
      dtalp -> (dtbetasq - dg4[0,0,0]) / (2 alp),

      k[li,lj] -> (- dtg[li,lj]
                   + g[lk,lj] dbeta[uk,li] + g[li,lk] dbeta[uk,lj]
                   + dg[li,lj,lk] beta[uk]) / (2 alp)
    } /. arrayToSymbolRules /. symmetryRules)
  };

  calculations =
  {
    calc["initial"],
    calc["always"]
  };

  CreateKrancThornTT[admGroups, "../", thorn,
    Calculations -> calculations,
    RealParameters -> Join[realParameters, parameters],
    ExtendedKeywordParameters -> extendedKeywordParameters,
    InheritedImplementations -> {"admbase"},
    CSE -> True,
    UseLoopControl -> True,
    (* We currently need to disable vectorisation, because it cannot
       handle the Piecewise[] functions that are generated from
       derivatives of the Min[] function used to hide singularities
       e.g. in ModifiedSchwarzschildBL *)
    UseVectors -> False];
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
thornNameRules = {"Vaidya" -> "Vaidya2"};

thorns = spacetimes /. thornNameRules;
MapThread[idThorn, {spacetimes, thorns}];

DeleteFile["../doc/spacetimes.tex"];
MapThread[docGenerate, {spacetimes, thorns}];
