#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 1: first-order singular mode, consistency with the 3D singular field (note sections 4, 6, 7).
   Run:  wolframscript -file stage1_singular_mode.wls                                        *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
SetCSParameters[<|"q" -> 1, "r0" -> 5|>];
r0 = csPar["r0"]; q = csPar["q"];
lmax = 6;            (* main range *)
lmaxPattern = 12;    (* extra range, used only to confirm the closed form of D_l *)
csResetChecks[];
t0 = AbsoluteTime[];

exactPts = {1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 49/10, 51/10, 11/2, 6, 13/2, 7, 8, 9, 10};
upperPts = Flatten@Table[x + I y, {x, {1, 24/5, 5, 26/5, 9}}, {y, {3/10, 1, 2}}];
realGrid = Select[Subdivide[1, 10, 180], Abs[# - r0] > 1/100 &];
aGrid = Subdivide[tubeL[], tubeR[], 100];
z = zVar;

Print["=== Stage 1: first-order singular mode (r0 = ", r0, ", contour amp = ", csPar["amp"], ", radius = ", csPar["radius"], ") ==="];

(* ---------------------------------------------------------------- *)
Print["\n--- 1a. What DSolve returns (psi_p = DSolveValue[L_l psi == deltaExt] /. C[_] -> 0) ---"];
Do[
  Print["l = ", l, ":  psi_p(z) = ", InputForm[Collect[Expand[dsolveSingularMode[l][z]], {Log[_]}, Simplify]]];
  Print["        logarithms: ", InputForm[Union[Cases[dsolveSingularMode[l][z], Log[_], Infinity]]]],
  {l, 0, lmax}];

(* ---------------------------------------------------------------- *)
Print["\n--- 1b. Branch handling ---"];
Do[
  Module[{raw = dsolveSingularMode[l], can = singularModeCanonical[l], dUHP, dAbove, dNaive},
   (* same analytic function in the upper half plane *)
   dUHP = Max[Abs[N[can /@ upperPts, 40] - N[raw /@ upperPts, 40]]];
   csCheck["l=" <> ToString[l] <> " canonical form == DSolve form in upper half plane", dUHP, 10^-30];
   (* principal-branch value of the canonical form on the axis == boundary value from above *)
   dAbove = Max[Abs[N[can /@ realGrid, 40] - (fromAbove[raw, #] & /@ realGrid)]];
   csCheck["l=" <> ToString[l] <> " canonical form on axis == raw form at r + i 1e-40", dAbove, 10^-30];
   (* pitfall: evaluating the raw DSolve form directly at real r *)
   dNaive = Max[Abs[N[raw /@ Select[realGrid, # > r0 &], 40] - (fromAbove[raw, #] & /@ Select[realGrid, # > r0 &])]];
   Print["   l=", l, ": raw DSolve form evaluated directly at real r > r0 differs from the from-above value by up to ", sci[dNaive],
    "  (Log[r0 - r] taken on the wrong side of its cut)"];
   ], {l, 0, lmax}];

(* ---------------------------------------------------------------- *)
Print["\n--- 1c. D_l = Re psi_p(r+i0) - G_l(r) ---"];
dRows = {};
Do[
  Module[{ab, al, be, res, abR, naive, naiveAB, naiveRes, note, mirror},
   ab = singularModeD[l]; {al, be} = ab;
   (* verify the exact two-parameter form at every exact point, both sides of r0 *)
   res = Max[Abs[N[(singularModeDValue[l, #] - (al #^l + be #^(-l - 1))) & /@ exactPts, 40]]];
   csCheck["l=" <> ToString[l] <> " D_l = alpha r^l + beta r^(-l-1) (fit at r=2,3; checked at 17 points, both sides)", res, 10^-30];
   (* independent fit on the right side only *)
   abR = Module[{d7 = singularModeDValue[l, 7], d8 = singularModeDValue[l, 8], A, B},
     {A, B} /. First@Solve[{A 7^l + B 7^(-l - 1) == d7, A 8^l + B 8^(-l - 1) == d8}, {A, B}]];
   csCheck["l=" <> ToString[l] <> " fit at r=7,8 gives the same (alpha,beta)", Max[Abs[N[abR - ab, 40]]], 10^-30];
   (* numerical cross-check with the raw form at r + i eps *)
   csCheck["l=" <> ToString[l] <> " D_l from raw form at r + i 1e-40 on a 180-point grid",
    Max[Abs[(Re[fromAbove[dsolveSingularMode[l], #]] - N[coulombModeUnit[l, #], 50] - N[al #^l + be #^(-l - 1), 50]) & /@ realGrid]], 10^-30];
   (* the naive (principal branch on the real axis) evaluation is NOT a homogeneous solution across r0 *)
   naive[r_] := Re[N[dsolveSingularMode[l][r], 40]] - N[coulombModeUnit[l, r], 40];
   naiveAB = Module[{A, B}, {A, B} /. First@Solve[{A 2^l + B 2^(-l - 1) == naive[2], A 3^l + B 3^(-l - 1) == naive[3]}, {A, B}]];
   naiveRes = Max[Abs[(naive[#] - (naiveAB[[1]] #^l + naiveAB[[2]] #^(-l - 1))) & /@ exactPts]];
   Print["   l=", l, ":  alpha = ", InputForm[al], ",  beta = ", InputForm[be],
    "    [naive real-axis evaluation: residual of the same fit on r > r0 = ", sci[naiveRes], "]"];
   note = {0, r0^(l + 2)/(2 l + 1)}; mirror = {r0^(1 - l)/(2 l + 1), 0};
   AppendTo[dRows, {l, al, be, note[[2]], mirror[[1]], Which[ab === note, "note (Log[z-r0])", ab === mirror, "mirror (Log[r0-z])", ab === note + mirror, "note + mirror", ab === {0, 0}, "D = 0", True, "other"]}];
   ], {l, 0, lmax}];

(* pattern check up to lmaxPattern *)
patternOK = And @@ Table[
    singularModeD[l] === Which[l == 0, {r0, r0^2}, l == 1, {0, 0}, True, {r0^(1 - l)/(2 l + 1), 0}], {l, 0, lmaxPattern}];
csCheck["closed form alpha_0=r0, beta_0=r0^2; D_1=0; alpha_l=r0^(1-l)/(2l+1), beta_l=0 (l>=2), for l<=" <> ToString[lmaxPattern],
  If[patternOK, 0, 1], 0];
Do[If[l > lmax, AppendTo[dRows, Module[{ab = singularModeD[l]}, {l, ab[[1]], ab[[2]], r0^(l + 2)/(2 l + 1), r0^(1 - l)/(2 l + 1),
      Which[ab === {r0^(1 - l)/(2 l + 1), 0}, "mirror (Log[r0-z])", True, "other"]}]]], {l, 0, lmaxPattern}];

Print["\n   Table: D_l = alpha r^l + beta r^(-l-1)  (note predicts alpha = 0, beta = r0^(l+2)/(2l+1); mirror choice alpha = r0^(1-l)/(2l+1), beta = 0)"];
Print["   l | alpha | beta | note beta | mirror alpha | matches"];
Scan[Print["   ", StringRiffle[ToString[#, InputForm] & /@ #, " | "]] &, dRows];
Export[FileNameJoin[{resDir, "stage1_D_table.csv"}],
  Prepend[Map[ToString[#, InputForm] &, dRows, {2}], {"l", "alpha", "beta", "note_beta_r0^(l+2)/(2l+1)", "mirror_alpha_r0^(1-l)/(2l+1)", "matches"}], "CSV"];

(* ---------------------------------------------------------------- *)
Print["\n--- 1d. Corrected singular mode phiS_l = psi_p - D_l: conditions (C1), (C2), continuity ---"];
c1Rows = {};
Do[
  Module[{phiS = singularModeCorrected[l], c1exact, c1num, c2sym, c2con, u, dCont, dEnd, dEndNaive, rawMinusD, imEnd},
   (* (C1) exact at rational points *)
   c1exact = Max[Abs[N[(ComplexExpand[Re[phiS[#]]] - coulombModeUnit[l, #]) & /@ exactPts, 40]]];
   csCheck["l=" <> ToString[l] <> " (C1) Re phiS(r+i0) = G_l exactly at 17 points", c1exact, 10^-30];
   (* (C1) numerically, from the raw DSolve form minus D at r + i eps *)
   rawMinusD = Function @@ {{zVar}, dsolveSingularMode[l][zVar] - singularModeD[l][[1]] zVar^l - singularModeD[l][[2]] zVar^(-l - 1)};
   c1num = Max[Abs[(Re[fromAbove[rawMinusD, #]] - N[coulombModeUnit[l, #], 50]) & /@ realGrid]];
   csCheck["l=" <> ToString[l] <> " (C1) on a 180-point grid via r + i 1e-40", c1num, 10^-30];
   (* (C2) symbolic *)
   c2sym = Simplify[radialL[l, phiS[z], z] - deltaExt[z]];
   csCheck["l=" <> ToString[l] <> " (C2) L_l phiS - deltaExt simplifies to 0", If[c2sym === 0, 0, 1], 0];
   (* (C2) numerically on the contour with the notebook's operator *)
   c2con = Module[{a, expr}, expr = contourL[l, Function @@ {{a}, phiS[yC[a]]}, a] - deltaExt[yC[a]];
     Max[Abs[(expr /. a -> N[#, 60]) & /@ Most[Rest[aGrid]]]]];
   csCheck["l=" <> ToString[l] <> " (C2) on the contour, odeOpComplex[phiS] - deltaExt", c2con, 10^-25];
   (* continuity along the contour: integrate L u = deltaExt along the contour from aL and compare *)
   u = solveContourSegment[l, Module[{a}, numericFunction[deltaExt[yC[a]], a]],
     {phiS[tubeL[]], D[phiS[z], z] /. z -> tubeL[]}];
   dCont = Max[Abs[(u[#] - N[phiS[yC[#]], 40]) & /@ aGrid]];
   csCheck["l=" <> ToString[l] <> " phiS(y(a)) equals the ODE continuation along the contour (101 points, interpolated)", dCont, 10^-12];
   dEnd = Abs[u[tubeR[]] - N[phiS[tubeR[]], 40]];
   csCheck["l=" <> ToString[l] <> " ODE continuation to aR equals phiS(aR + i0) (canonical form)", dEnd, 10^-15];
   (* the same comparison at aR using the raw DSolve form evaluated with Mathematica's principal branch *)
   dEndNaive = Abs[u[tubeR[]] - N[rawMinusD[tubeR[]], 40]];
   imEnd = Im[u[tubeR[]]];
   Print["   l=", l, ": at aR, ODE continuation vs canonical phiS(aR) = ", sci[dEnd],
    ";  vs raw form evaluated at real aR (principal branch) = ", sci[dEndNaive]];
   AppendTo[c1Rows, {l, c1exact, c1num, c2sym === 0, c2con, dCont, dEndNaive}];
   ], {l, 0, lmax}];
Export[FileNameJoin[{resDir, "stage1_conditions.csv"}],
  Prepend[Map[If[NumericQ[#], sci[#], ToString[#]] &, c1Rows, {2}],
   {"l", "C1_exact_max", "C1_numeric_max", "C2_symbolic_zero", "C2_contour_max", "contour_continuity_max", "aR_raw_principal_branch_error"}], "CSV"];

(* ---------------------------------------------------------------- *)
Print["\n--- 1e. 3D meaning of the uncorrected psi_p ---"];
(* Delta(x) = Sum_l kappa_l D_l(r) Y_l(theta).  With the closed form above, for r < r0:
   Delta = -(q/r0) Sum_{l != 1} (r/r0)^l P_l - q/r = -q/rho + q r cos(theta)/r0^2 - q/r ; divergent for r > r0. *)
Module[{rr = 2, th = 1, L = 80, sum, closed},
  sum = Sum[kappaL[l] (Which[l == 0, r0 + r0^2/rr, l == 1, 0, True, r0^(1 - l) rr^l/(2 l + 1)]) Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, Cos[th]], {l, 0, L}];
  closed = -q/Sqrt[rr^2 + r0^2 - 2 rr r0 Cos[th]] + q rr Cos[th]/r0^2 - q/rr;
  csCheck["3D sum of kappa_l D_l Y_l at r=2, theta=1 (l<=80) equals -q/rho + q r cos(theta)/r0^2 - q/r", N[sum - closed, 40], 10^-25]];

(* ---------------------------------------------------------------- *)
Print["\n--- 1f. Plots ---"];
plots = Table[
   Module[{raw = dsolveSingularMode[l], phiS = singularModeCorrected[l]},
    Plot[{Re[N[singularModeCanonical[l][r]]], Re[N[raw[r]]], coulombModeUnit[l, r], Re[N[phiS[r]]]}, {r, 1, 10},
     PlotStyle -> {{Blue, Thick}, {Orange, Dashed}, {Black, Thick}, {Red, Dotted, Thick}},
     PlotLegends -> {"Re \[Psi]p(r+i0)", "Re \[Psi]p naive (principal branch at real r)", "\!\(\*SubscriptBox[\(G\), \(\[ScriptL]\)]\)", "Re \!\(\*SuperscriptBox[\(\[Phi]\), \(S\)]\)(r+i0) corrected"},
     PlotLabel -> "\[ScriptL] = " <> ToString[l], AxesLabel -> {"r", ""}, ImageSize -> 420, PlotRange -> All, Exclusions -> {r == r0}]],
   {l, 0, 2}];
Export[FileNameJoin[{resDir, "stage1_singular_modes.pdf"}], Column[plots]];
contourPlots = Table[
   Module[{phiS = singularModeCorrected[l]},
    Plot[{Re[N[phiS[yC[a]]]], Im[N[phiS[yC[a]]]]}, {a, tubeL[], tubeR[]}, PlotStyle -> {Blue, Red},
     PlotLegends -> {"Re \!\(\*SuperscriptBox[\(\[Phi]\), \(S\)]\)(y(a))", "Im \!\(\*SuperscriptBox[\(\[Phi]\), \(S\)]\)(y(a))"},
     PlotLabel -> "corrected singular mode on the contour, \[ScriptL] = " <> ToString[l], AxesLabel -> {"a", ""}, ImageSize -> 420]],
   {l, 0, 2}];
Export[FileNameJoin[{resDir, "stage1_phiS_on_contour.pdf"}], Column[contourPlots]];
Print["   wrote stage1_singular_modes.pdf, stage1_phiS_on_contour.pdf, stage1_D_table.csv, stage1_conditions.csv"];

Print["\nStage 1 wall time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
allPass = csSummary[];
Export[FileNameJoin[{resDir, "stage1_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
