#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 5a: the 1/rho^4 case, Q[A,B] = Grad A . Grad B.  Source modes, couplings, continuations and analytic modes.
   PhiR = c0 + c1 r cos(theta) + c2 r^2 P2(cos(theta)), c0 = 1, c1 = 1/10, c2 = 1/100.
   Run:  wolframscript -file stage5a_gradient_sources.wls                                        *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10, "c2" -> 1/100|>];
Q = "GradDot";
r0 = csPar["r0"]; q = csPar["q"]; c0 = csPar["c0"]; c1 = csPar["c1"]; c2 = csPar["c2"]; aL = tubeL[]; aR = tubeR[];
lmax = 6;
csResetChecks[];
t0 = AbsoluteTime[];
testR = {1, 2, 4, 9/2, 49/10, 51/10, 11/2, 7, 10};
nint[f_] := Quiet@NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> 50, PrecisionGoal -> 25, AccuracyGoal -> 35, MaxRecursion -> 30];  (* AccuracyGoal: some integrals vanish *)
Y[l_, mu_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, mu];
relDiff[a_, b_] := Abs[a - b]/Max[Abs[b], 10^-30];
rhoF[r_, mu_] := Sqrt[r^2 + r0^2 - 2 r r0 mu];

Print["=== Stage 5a: 1/rho^4 case (Q = Grad.Grad): sources, couplings, analytic modes ==="];

(* ---------------- SS ---------------- *)
Print["\n--- SS: [1/rho^4]_l = -Sqrt[Pi(2l+1)] Q_l'(x)/(2 r^2 r0^2) ---"];
csCheck["l=0 closed form 2 Sqrt[Pi]/(r^2 - r0^2)^2 at 9 radii",
  Max[Table[relDiff[N[ssSourceModeExpr[Q, 0, r, "real"], 40], N[q^2 2 Sqrt[Pi]/(r^2 - r0^2)^2, 40]], {r, testR}]], 10^-30];
Do[
  csCheck["l=" <> ToString[l] <> " Q_l'-formula vs direct mu-integration of 1/rho^4 (rel)",
   Max[Table[relDiff[N[ssSourceModeExpr[Q, l, r, "real"], 40], q^2 Sqrt[Pi (2 l + 1)] nint[LegendreP[l, mu]/rhoF[r, mu]^4]], {r, testR}]], 10^-20];
  csCheck["l=" <> ToString[l] <> " Q_l' (Neumann form) equals D[LegendreQ[l,0,3,x], x]",
   Max[Table[Module[{x = xOfZ[r]}, relDiff[N[legendreQNeumannDx[l, x, Log[(x + 1)/(x - 1)], 1/(x^2 - 1)], 40],
       N[D[LegendreQ[l, 0, 3, xx], xx] /. xx -> N[x, 50], 40]]], {r, testR}]], 10^-25];
  csCheck["l=" <> ToString[l] <> " Q_l'-formula vs closed-form u-integral (rel)",
   Max[Table[relDiff[N[ssSourceModeExpr[Q, l, r, "real"], 40], analyticEval[pieceAnalyticSource[Q, "SS", l], r]], {r, testR}]], 10^-25];
  csCheck["l=" <> ToString[l] <> " Re SS_continued(r+i0) = SS_real(r) for r != r0 (rel)",
   Max[Table[relDiff[Re[N[ssSourceModeExpr[Q, l, r, "continued"], 40]], N[ssSourceModeExpr[Q, l, r, "real"], 40]], {r, testR}]], 10^-28];
  Module[{f = ssSourceContour[Q, l], vals},
   vals = f /@ N[Subdivide[aL, aR, 4000], 40];
   csCheck["l=" <> ToString[l] <> " SS continued source along the contour: max |step|/max|f| (no jumps)", Max[Abs[Differences[vals]]]/Max[Abs[vals]], 10^-2]];
  (* the identity behind the regularisation test: continued source = L_z of the continued (q^2/2)[1/rho^2]_l, identically *)
  Module[{z, psi},
   psi = (1/2) ssSourceModeExpr["Product", l, z, "continued"];      (* (q^2/2) [1/rho^2]_l, continued with Log[z - r0] *)
   csCheck["l=" <> ToString[l] <> " L_z[(q^2/2)[1/rho^2]_l continued] - [q^2/rho^4]_l continued == 0 (symbolic)",
    If[Simplify[radialL[l, psi, z] - ssSourceModeExpr[Q, l, z, "continued"]] === 0, 0, 1], 0]],
  {l, 0, lmax}];

(* ---------------- couplings ---------------- *)
Print["\n--- SR couplings d_r K and K^theta ---"];
Module[{pairs = {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0, 2}, {2, 0}, {1, 2}, {2, 1}, {2, 2}, {3, 2}, {4, 2}, {5, 1}, {6, 2}}, dK, dT},
  dK = Max[Flatten@Table[relDiff[N[KdCoupling[p[[1]], p[[2]], r], 40],
       2 Pi nint[Y[p[[1]], mu] Y[p[[2]], mu] (-(r - r0 mu)/rhoF[r, mu]^3)]], {p, pairs}, {r, {2, 47/10, 53/10, 8}}]];
  csCheck["d_r K_{l l'} vs direct angular integration of d_r(1/rho) Y_l Y_l' (13 pairs x 4 radii, rel)", dK, 10^-20];
  dT = Max[Flatten@Table[relDiff[N[KthetaCoupling[p[[1]], p[[2]], r], 40],
       2 Pi nint[(1 - mu^2) (r r0/rhoF[r, mu]^3) D[Y[p[[2]], mu], mu] Y[p[[1]], mu]]], {p, pairs}, {r, {2, 47/10, 53/10, 8}}]];
  csCheck["K^theta_{l l'} = Sum_L <l l' L>[1/rho]_L (L(L+1)+l'(l'+1)-l(l+1))/2 vs direct integration of d_theta(1/rho) d_theta Y_l' Y_l (rel)", dT, 10^-20]];
Do[
  csCheck["l=" <> ToString[l] <> " Re of continued d_r K, K^theta (r+i0) = real values, Theta->H and corrected-phiS",
   Max[Flatten@Table[{Abs[Re[N[KdCouplingContinued[l, lp, r, m], 40]] - N[KdCoupling[l, lp, r], 40]],
       Abs[Re[N[KthetaCouplingContinued[l, lp, r, m], 40]] - N[KthetaCoupling[l, lp, r], 40]]}, {lp, 0, 2}, {r, DeleteCases[testR, r0]}, {m, {"ThetaToH", "CorrectedSingular"}}]], 10^-28],
  {l, 0, lmax}];

(* ---------------- SR, RR source modes ---------------- *)
Print["\n--- SR and RR source modes ---"];
phiRpair = Association@Table[With[{f = phiRAnalytic[lp]}, lp -> {f, f'}], {lp, 0, 2}];
Do[
  Module[{sr = srSourceReal[Q, l, phiRpair], rr = rrSourceReal[Q, l, phiRpair], solSR = pieceAnalyticSolution[Q, "SR", l], pp, J, p0 = phiRAtParticle[]},
   csCheck["l=" <> ToString[l] <> " SR coupling formula vs closed form of [2q Grad(1/rho).Grad PhiR]_l (rel)",
    Max[Table[relDiff[N[sr[r], 40], analyticEval[pieceAnalyticSource[Q, "SR", l], r]], {r, DeleteCases[testR, r0]}]], 10^-25];
   csCheck["l=" <> ToString[l] <> " SR closed form vs direct mu-integration of the 3D source (rel)",
    Max[Table[relDiff[analyticEval[pieceAnalyticSource[Q, "SR", l], r],
       Sqrt[Pi (2 l + 1)] nint[LegendreP[l, mu] (-2 q (-c2 r^2 (1 - mu^2) + (r mu - r0) (c1 + 2 c2 r mu)))/rhoF[r, mu]^3]], {r, testR}]], 10^-18];
   (* identity 2 Grad PhiS.Grad PhiR = Lap(PhiS PhiR) + 4 Pi q PhiR(x0) delta^3: pointwise part and the delta(r - r0) part *)
   pp = Table[uModeIntegral[l, 0, "invSqrt", zVar, br] q (c0 - c2 zVar^2/2) + q c1 zVar uModeIntegral[l, 1, "invSqrt", zVar, br] +
       q (3/2) c2 zVar^2 uModeIntegral[l, 2, "invSqrt", zVar, br], {br, {"in", "out"}}];      (* [q PhiR/rho]_l *)
   csCheck["l=" <> ToString[l] <> " identity, pointwise: L_l [q PhiR/rho]_l = SR source mode for r != r0 (symbolic, both branches)",
    If[Simplify[{radialL[l, pp[[1]], zVar], radialL[l, pp[[2]], zVar]} - pieceAnalyticSource[Q, "SR", l], zVar > 0] === {0, 0}, 0, 1], 0];
   J = Limit[D[pp[[2]], zVar], zVar -> r0] - Limit[D[pp[[1]], zVar], zVar -> r0];
   csCheck["l=" <> ToString[l] <> " identity, delta part: jump of d_r [q PhiR/rho]_l at r0 = -4 Pi q PhiR(x0) Y_l(0)/r0^2",
    N[Simplify[J + 4 Pi q p0 Sqrt[(2 l + 1)/(4 Pi)]/r0^2], 30], 10^-25];
   csCheck["l=" <> ToString[l] <> " RR gradient-Gaunt coupling vs closed form of [|Grad PhiR|^2]_l",
    Max[Table[Abs[N[rr[r], 40] - analyticEval[pieceAnalyticSource[Q, "RR", l], r]], {r, testR}]], 10^-25]],
  {l, 0, lmax}];

(* ---------------- analytic particular solutions ---------------- *)
Print["\n--- Analytic particular solutions: q^2/(2 rho^2), q (PhiR - PhiR(x0))/rho, PhiR^2/2 ---"];
Do[
  Do[
   Module[{sol = pieceAnalyticSolution[Q, piece, l], src = pieceAnalyticSource[Q, piece, l], res, cont, direct},
    res = Simplify[{radialL[l, sol[[1]], zVar] - src[[1]], radialL[l, sol[[2]], zVar] - src[[2]]}, zVar > 0];
    csCheck[piece <> " l=" <> ToString[l] <> " L_l(solution mode) - source mode == 0 (symbolic, both branches)", If[res === {0, 0}, 0, 1], 0];
    If[piece =!= "SS",
     cont = Abs[Limit[sol[[1]], zVar -> r0] - Limit[sol[[2]], zVar -> r0]] + Abs[Limit[D[sol[[1]], zVar], zVar -> r0] - Limit[D[sol[[2]], zVar], zVar -> r0]];
     csCheck[piece <> " l=" <> ToString[l] <> " solution mode is C^1 across r0", N[cont, 30], 10^-25]];
    direct = Max[Table[relDiff[analyticEval[sol, r], Sqrt[Pi (2 l + 1)] nint[analytic3D[Q, piece, r, ArcCos[mu]] LegendreP[l, mu]]], {r, {2, 47/10, 53/10, 8}}]];
    csCheck[piece <> " l=" <> ToString[l] <> " solution mode vs direct mu-integration of the 3D function (rel)", direct, 10^-18]],
   {piece, {"SS", "SR", "RR"}}], {l, 0, lmax}];

Print["\nStage 5a wall time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage5a_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
