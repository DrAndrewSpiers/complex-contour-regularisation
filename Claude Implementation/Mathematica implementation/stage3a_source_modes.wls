#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 3a: second-order source modes, couplings, continuations and analytic particular solutions (1/rho^2 case).
   Run:  wolframscript -file stage3a_source_modes.wls                                        *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>];
r0 = csPar["r0"]; q = csPar["q"]; c0 = csPar["c0"]; c1 = csPar["c1"]; aL = tubeL[]; aR = tubeR[];
lmax = 6;
csResetChecks[];
t0 = AbsoluteTime[];
z = zVar;
testR = {1, 2, 4, 9/2, 49/10, 51/10, 11/2, 7, 10};          (* both sides of r0, inside and outside the tube *)
nint[f_] := Quiet@NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 30];  (* Quiet: slow-convergence warnings for integrals that vanish *)
Y[l_, mu_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, mu];
relDiff[a_, b_] := Abs[a - b]/Max[Abs[b], 10^-30];

Print["=== Stage 3a: source modes and couplings (Q = Product, 1/rho^2 case) ==="];

(* ---------------------------------------------------------------- *)
Print["\n--- SS: Neumann formula [1/rho^2]_l = Sqrt[Pi(2l+1)] Q_l(x)/(r r0), x = (r^2+r0^2)/(2 r r0) ---"];
Do[
  Module[{dN, dQ3, dU},
   dN = Max[Table[relDiff[N[ssSourceModeExpr["Product", l, r, "real"], 40],
       Sqrt[Pi (2 l + 1)] nint[LegendreP[l, mu]/(r^2 + r0^2 - 2 r r0 mu)]], {r, testR}]];
   dQ3 = Max[Table[relDiff[N[legendreQNeumann[l, xOfZ[r], Log[(xOfZ[r] + 1)/(xOfZ[r] - 1)]], 40], LegendreQ[l, 0, 3, N[xOfZ[r], 40]]], {r, testR}]];
   dU = Max[Table[relDiff[N[ssSourceModeExpr["Product", l, r, "real"], 40], analyticEval[pieceAnalyticSource["Product", "SS", l], r]], {r, testR}]];
   csCheck["l=" <> ToString[l] <> " Neumann form vs direct mu-integration (rel)", dN, 10^-20];
   csCheck["l=" <> ToString[l] <> " Q_l via P_l log/2 - W_(l-1) equals LegendreQ[l,0,3,x] (type 3)", dQ3, 10^-25];
   csCheck["l=" <> ToString[l] <> " Neumann form vs closed-form u-integral (rel)", dU, 10^-25]], {l, 0, lmax}];

Print["\n--- SS: continuation 2 Log[z+r0] - 2 Log[z-r0]: Re on the axis (from above) = real source ---"];
Do[
  Module[{dRe, imL, imR, jump, fA, vals},
   dRe = Max[Table[relDiff[Re[N[ssSourceModeExpr["Product", l, r, "continued"], 40]], N[ssSourceModeExpr["Product", l, r, "real"], 40]], {r, testR}]];
   csCheck["l=" <> ToString[l] <> " Re SS_continued(r+i0) = SS_real(r) at 9 points (rel)", dRe, 10^-28];
   (* continuity along the contour: no 2 Pi i jumps *)
   fA = ssSourceContour["Product", l];
   vals = fA /@ N[Subdivide[aL, aR, 4000], 40];
   jump = Max[Abs[Differences[vals]]]/Max[Abs[vals]];
   csCheck["l=" <> ToString[l] <> " SS continued source along contour: max |step|/max|f| on 4000 steps (no jumps)", jump, 10^-2];
   imL = Im[N[ssSourceModeExpr["Product", l, 47/10, "continued"], 20]]; imR = Im[N[ssSourceModeExpr["Product", l, 53/10, "continued"], 20]];
   If[l <= 2, Print["   l=", l, ": shadow Im SS(r+i0) at r=4.7: ", sci[imL], ", at r=5.3: ", sci[imR], "  (Im = -Pi P_l(x) prefactor for r<r0, 0 for r>r0)"]]], {l, 0, lmax}];

(* ---------------------------------------------------------------- *)
Print["\n--- Couplings: Gaunt coefficients and K_{l l'} ---"];
csCheck["Gaunt <l1 l2 l3> vs exact integration 2 Pi Integrate[Y Y Y], all l1,l2 <= 4, l3 <= 8",
  Max[Flatten@Table[Abs[Simplify[gaunt[l1, l2, l3] - 2 Pi Integrate[Y[l1, mu] Y[l2, mu] Y[l3, mu], {mu, -1, 1}]]], {l1, 0, 4}, {l2, 0, 4}, {l3, 0, 8}]], 0];
Module[{pairs = {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {2, 1}, {3, 1}, {4, 0}, {5, 1}, {6, 1}, {2, 2}}, d},
  d = Max[Flatten@Table[relDiff[N[KCoupling[p[[1]], p[[2]], r], 40], 2 Pi nint[Y[p[[1]], mu] Y[p[[2]], mu]/Sqrt[r^2 + r0^2 - 2 r r0 mu]]], {p, pairs}, {r, {2, 47/10, 53/10, 8}}]];
  csCheck["K_{l l'}(r) vs direct angular integration (10 pairs x 4 radii, rel)", d, 10^-20]];
Do[
  Module[{dA, dB, imA, imB},
   dA = Max[Flatten@Table[Abs[Re[N[KCouplingContinued[l, lp, r, "ThetaToH"], 40]] - N[KCoupling[l, lp, r], 40]], {lp, 0, 1}, {r, testR}]];
   dB = Max[Flatten@Table[Abs[Re[N[KCouplingContinued[l, lp, r, "CorrectedSingular"], 40]] - N[KCoupling[l, lp, r], 40]], {lp, 0, 1}, {r, testR}]];
   csCheck["l=" <> ToString[l] <> " Re K_continued(r+i0) = K(r): Theta->H", dA, 10^-28];
   csCheck["l=" <> ToString[l] <> " Re K_continued(r+i0) = K(r): corrected singular mode (Stage 1)", dB, 10^-28];
   imA = Im[N[KCouplingContinued[l, 1, 47/10, "ThetaToH"], 20]]; imB = Im[N[KCouplingContinued[l, 1, 47/10, "CorrectedSingular"], 20]];
   If[l <= 2, Print["   l=", l, ", l'=1, r=4.7: Im K (Theta->H) = ", sci[imA], ", Im K (corrected) = ", sci[imB], "  (shadows differ, as expected)"]]], {l, 0, lmax}];

(* ---------------------------------------------------------------- *)
Print["\n--- SR and RR source modes: coupling formula vs closed-form angular integrals ---"];
phiRexact = <|0 -> phiRAnalytic[0], 1 -> phiRAnalytic[1]|>;
Do[
  Module[{srK, rrG},
   srK = srSourceReal["Product", l, phiRexact];
   csCheck["l=" <> ToString[l] <> " SR: 2q Sum K_{l l'} phiR_l' vs [2 q PhiR/rho]_l closed form (rel)",
    Max[Table[relDiff[N[srK[r], 40], analyticEval[pieceAnalyticSource["Product", "SR", l], r]], {r, testR}]], 10^-25];
   csCheck["l=" <> ToString[l] <> " SR closed form vs direct mu-integration (rel)",
    Max[Table[relDiff[analyticEval[pieceAnalyticSource["Product", "SR", l], r],
       Sqrt[Pi (2 l + 1)] nint[2 q (c0 + c1 r mu) LegendreP[l, mu]/Sqrt[r^2 + r0^2 - 2 r r0 mu]]], {r, testR}]], 10^-20];
   rrG = rrSourceReal["Product", l, phiRexact];
   csCheck["l=" <> ToString[l] <> " RR: Gaunt coupling vs [(PhiR)^2]_l closed form",
    Max[Table[Abs[N[rrG[r], 40] - analyticEval[pieceAnalyticSource["Product", "RR", l], r]], {r, testR}]], 10^-25]], {l, 0, lmax}];

(* ---------------------------------------------------------------- *)
Print["\n--- Analytic particular solutions (modes of q^2 Log rho, the SR check, the RR check) ---"];
Module[{greens = (Sqrt[Pi] (-4 zVar r0 - (zVar - r0)^2 Log[(zVar - r0)^2] + (zVar + r0)^2 Log[(zVar + r0)^2]))/(4 zVar r0)},
  csCheck["SS l=0 mode of Log rho equals AnalyticSol of Greens-function-Poisson-equation-3.nb",
   Max[Table[Abs[analyticEval[pieceAnalyticSolution["Product", "SS", 0], r] - N[greens /. zVar -> r, 30]], {r, testR}]], 10^-25]];
Do[
  Do[
   Module[{sol = pieceAnalyticSolution["Product", piece, l], src = pieceAnalyticSource["Product", piece, l], res, cont, direct},
    res = Simplify[{radialL[l, sol[[1]], zVar] - src[[1]], radialL[l, sol[[2]], zVar] - src[[2]]}, zVar > 0];
    csCheck[piece <> " l=" <> ToString[l] <> " L_l(solution mode) - source mode == 0 (symbolic, both branches)", If[res === {0, 0}, 0, 1], 0];
    (* C^1 matching across r0 of the two branch formulas *)
    cont = Abs[Limit[sol[[1]], zVar -> r0] - Limit[sol[[2]], zVar -> r0]] + Abs[Limit[D[sol[[1]], zVar], zVar -> r0] - Limit[D[sol[[2]], zVar], zVar -> r0]];
    csCheck[piece <> " l=" <> ToString[l] <> " solution mode is C^1 across r0 (no point-supported source)", N[cont, 30], 10^-25];
    direct = Max[Table[relDiff[analyticEval[sol, r],
        Sqrt[Pi (2 l + 1)] nint[analytic3D[piece, r, ArcCos[mu]] LegendreP[l, mu]]], {r, {2, 47/10, 53/10, 8}}]];
    csCheck[piece <> " l=" <> ToString[l] <> " solution mode vs direct mu-integration of the 3D function (rel)", direct, 10^-18]],
   {piece, {"SS", "SR", "RR"}}], {l, 0, lmax}];

Print["\n(The Grad.Grad case, Q = \"GradDot\", is tested in stage5a/stage5b.)"];

Print["\nStage 3a wall time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage3a_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
