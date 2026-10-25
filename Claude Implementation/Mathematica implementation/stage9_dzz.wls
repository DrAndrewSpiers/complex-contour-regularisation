#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 9: distributional SR.  Q[A,B] = (1/2) d_z^2 (A B), so Lap Phi2 = (1/2) d_z^2 (Phi1^2); PhiR = c0 + c1 z.
   Every Dzz piece is (1/2) d_z^2 of the Product piece; exact answer Phi2 = (1/2) d_z^2 chi (chi = Stage 3 total).
   The SR modes contain delta(r - r0) terms (d_z^2 of the kinked modes of 1/rho).
   Test 1: exact puncture q/rho; three-region solves (SS: regularisation; SR: product route, delta via the continuation).
   Test 2: Stage 6 puncture (n = 3, 5; R = 1/2): product route, H0 (Cauchy continuation of the pointwise real-line SR),
           H1(k), k = 0..3 (Part 1 = 2Q[PhiP, T_k] closed form continued; Part 2 = pointwise real line + Cauchy).
           Lemma-1 decomposition on the real axis, validated by contour solves.
   Run:  wolframscript -file stage9_dzz.wls [smoke]                                                                     *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
s8 = Import[FileNameJoin[{dir, "stage8_hybrid_split.wls"}], "Text"];
s8 = StringTake[s8, {1, StringPosition[s8, "(* ================= main loop ================= *)"][[1, 1]] - 1}];
s8 = StringReplace[s8, "dir = DirectoryName[$InputFileName];" -> "", 1];
s8 = StringReplace[s8, "Print[\"=== Stage 8" -> "Hold[\"x", 1];
ToExpression[s8];                                         (* Stage 6 + 7 + 8 definitions (smoke flag from Stage 7) *)
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
SetCSParameters[casePars["Product"]]; aL = tubeL[]; aR = tubeR[];          (* PhiR = c0 + c1 z, R = 1/2 *)
q = csPar["q"]; c0 = csPar["c0"]; c1 = csPar["c1"];
lOut9 = Range[0, 6]; lpList9 = {4, 8, 12, 16, 20, 30, 40}; nS9 = {3, 5}; kList9 = {0, 1, 2, 3}; Mprod = 20;
If[smoke, lOut9 = {0, 2}; lpList9 = {4, 8}; nS9 = {3}];
jMax = Max[lOut9] + 2; Mm = Max[lpList9]; lpProd = Select[lpList9, # <= Mprod &];
Print["=== Stage 9: Q = (1/2) d_z^2 (A B); l = ", lOut9, ", l'max = ", lpList9, ", n = ", nS9, " ==="];

(* response to delta(r - s) with zero data below s, and its r-derivative *)
G0[l_, r_, s_] := If[r > s, (s^2/(2 l + 1)) (r^l s^(-l - 1) - s^l r^(-l - 1)), 0];
dG0[l_, r_, s_] := If[r > s, (s^2/(2 l + 1)) (l r^(l - 1) s^(-l - 1) + (l + 1) s^l r^(-l - 2)), 0];
(* derivative lists: F_j -> {F, F', F'', ...}; Leibniz product; d_z on lists (output valid for j <= jtop - 1) *)
prodD[a_, b_] := Table[Sum[Binomial[m, i] a[[i + 1]] b[[m - i + 1]], {i, 0, m}], {m, 0, Length[a] - 1}];
invZD[p_, z_] := (-1)^p p!/z^(p + 1);
dzArr[Fl_, z_] := Module[{nd = Length[Fl[[1]]] - 1, f},
   f[j_, m_] := If[j < 0 || j > Length[Fl] - 1, 0, Fl[[j + 1, m + 1]]];
   Table[Table[aZ[j] (f[j - 1, m + 1] - (j - 1) Sum[Binomial[m, i] f[j - 1, i] invZD[m - i, z], {i, 0, m}]) +
      aZ[j + 1] (f[j + 1, m + 1] + (j + 2) Sum[Binomial[m, i] f[j + 1, i] invZD[m - i, z], {i, 0, m}]), {m, 0, nd - 1}], {j, 0, Length[Fl] - 2}]];
dzzArr[Fl_, z_] := dzArr[dzArr[Fl, z], z];
(* delta(r - r0) coefficient of (d_z^2 F)_l from one-sided first derivatives (F continuous) *)
dzzDeltaNum[l_, dIn_, dOut_] := Module[{J},
   J[j_] := If[j < 0, 0, If[j >= 1, aZ[j] (dOut[j - 1] - dIn[j - 1]), 0] + aZ[j + 1] (dOut[j + 1] - dIn[j + 1])];
   If[l >= 1, aZ[l] J[l - 1], 0] + aZ[l + 1] J[l + 1]];
invR[j_, R_] := Sqrt[4 Pi/(2 j + 1)] {zVar^j/R^(j + 1), R^j/zVar^(j + 1)};                 (* [1/rho]_j in/out, rho = |x - R z^| *)

(* ================= checks of the d_z machinery ================= *)
Off[N::meprec];         (* exact zeros in the symbolic checks *)
Do[Module[{f = Table[With[{b = b}, dzExpr[Function[j, If[j < 0, 0, invR[j, r0][[b]]]], l]], {b, 2}],
     ex = Table[-(zVar uModeIntegral[l, 1, "invSqrt3", zVar, b] - r0 uModeIntegral[l, 0, "invSqrt3", zVar, b]), {b, {"in", "out"}}], dr0, arr},
    csCheck["l=" <> ToString[l] <> ": d_z mode formula on 1/rho = modes of -(z - z0)/rho^3 (r = 3, 7)", Abs[N[f[[1]] - ex[[1]] /. zVar -> 3, 30]] + Abs[N[f[[2]] - ex[[2]] /. zVar -> 7, 30]], 10^-25];
    (* d_z^2 (1/rho) = d_r0^2 (1/rho): its delta coefficient is (d_r0 in - d_r0 out) at r = r0 *)
    dr0 = Module[{R0}, (D[invR[l, R0], R0] /. R0 -> r0) /. zVar -> r0];
    csCheck["l=" <> ToString[l] <> ": delta coefficient of d_z^2 [1/rho]_l = [d_r0 (in - out)](r0) = -4 Pi Y_l(0)/r0^2",
     Abs[N[dzzDeltaCoeff[invR[#, r0] &, l] - dr0 . {1, -1}, 30]] + Abs[N[dr0 . {1, -1} + 4 Pi Y0[l]/r0^2, 30]], 10^-25];
    (* numeric list d_z^2 against symbolic dzzExpr *)
    arr = Table[N[{#, D[#, zVar], D[#, {zVar, 2}], D[#, {zVar, 3}]} &[invR[j, r0][[2]]] /. zVar -> 7, 40], {j, 0, l + 2}];
    csCheck["l=" <> ToString[l] <> ": list d_z^2 (value, derivative) = symbolic dzzExpr at r = 7",
     Max[Abs[dzzArr[arr, 7][[l + 1, 1 ;; 2]] - N[{#, D[#, zVar]} &[dzzExpr[Function[j, If[j < 0, 0, invR[j, r0][[2]]]], l]] /. zVar -> 7, 40]]], 10^-30]], {l, 0, Max[lOut9]}];
Do[Do[Module[{sol = pieceAnalyticSolution["Dzz", p, l], src = pieceAnalyticSource["Dzz", p, l]},
     csCheck[p <> " l=" <> ToString[l] <> ": L_l((1/2) d_z^2 sol) = (1/2) d_z^2 src at r = 2, 7/2, 7, 9", Max[Table[Abs[N[(radialL[l, sol[[If[r < r0, 1, 2]]], zVar] - src[[If[r < r0, 1, 2]]]) /. zVar -> r, 30]], {r, {2, 7/2, 7, 9}}]], 10^-20]],
   {p, {"SS", "SR", "RR"}}], {l, lOut9}];
On[N::meprec];
Print["    machinery checks done [", tick[], " s]"];

(* ================= TEST 1: exact puncture q/rho ================= *)
Print["\n--- Test 1: exact puncture q/rho ---"];
(* product route on the contour: (1/2) d_z^2 of 2 q K~_{j l'} phiR_l' (K~: Theta -> H continuation; phiR analytic, entire) *)
srContExpr[l_] := (1/2) dzzExpr[Function[j, If[j < 0, 0, 2 q Sum[KCouplingContinued[j, lp, zVar, "ThetaToH"] phiRAnalytic[lp][zVar], {lp, 0, 1}]]], l];
srPwExpr[l_, b_] := (1/2) dzzExpr[Function[j, If[j < 0, 0, 2 q Sum[Sum[gaunt[j, lp, L] invR[L, r0][[b]], {L, Abs[j - lp], j + lp}] phiRAnalytic[lp][zVar], {lp, 0, 1}]]], l];
fitPoint9[l_, d_] := Module[{b1, b2, M},
   b1 = Function @@ {{zVar}, invR[l, r0][[2]] - invR[l, r0][[1]]};                          (* response to a point charge at x0 *)
   b2 = Function @@ {{zVar}, Module[{R0}, D[invR[l, R0][[2]] - invR[l, R0][[1]], R0] /. R0 -> r0]};   (* ... to a dipole *)
   M = N[{{b1[aR], b2[aR]}, {b1'[aR], b2'[aR]}}, 40]; LinearSolve[M, d]];
t1Rows = {};
Do[Module[{solSS, solSR, solRR, chiSol, scale, rS, rR, lR, lPw, fi, fo, srcR, eSS, eSR, eTot, ePw, eFix, abS, abR, gam},
    {solSS, solSR, solRR} = pieceAnalyticSolution["Dzz", #, l] & /@ {"SS", "SR", "RR"};
    chiSol = solSS + solSR + solRR;
    scale = Max[Abs[analyticEval[chiSol, #] & /@ gCom]];
    csCheck["Test 1 l=" <> ToString[l] <> ": product form of the pointwise SR = analytic (1/2) d_z^2 [2 q PhiR/rho] (r = 2, 7)",
     Max[Table[Abs[N[(srPwExpr[l, b] - pieceAnalyticSource["Dzz", "SR", l][[b]]) /. zVar -> {2, 7}[[b]], 30]], {b, 2}]], 10^-25];
    rS = threeRegionSolve[l, numericFunction[ssSourceModeExpr["Dzz", l, zVar, "real"], zVar], contourFunction[ssSourceModeExpr["Dzz", l, zVar, "continued"], zVar],
      {analyticEval[solSS, 1], analyticEvalD[solSS, 1]}];
    fi = numericFunction[pieceAnalyticSource["Dzz", "SR", l][[1]], zVar]; fo = numericFunction[pieceAnalyticSource["Dzz", "SR", l][[2]], zVar];
    srcR[r_?NumericQ] := If[r < r0, fi[r], fo[r]];
    rR = threeRegionSolve[l, srcR, contourFunction[srContExpr[l], zVar], {analyticEval[solSR, 1], analyticEvalD[solSR, 1]}];
    lR = realLineSolve[l, numericFunction[pieceAnalyticSource["Dzz", "RR", l][[1]], zVar], {analyticEval[solRR, 1], analyticEvalD[solRR, 1]}];
    lPw = lineSolve[l, fi, fo, {analyticEval[solSR, 1], analyticEvalD[solSR, 1]}];            (* negative control: pointwise SR, no delta *)
    gam = (1/2) N[dzzDeltaCoeff[pieceAnalyticSource["Product", "SR", #] &, l], 40];
    eSS = (exteriorValue[rS, #] - analyticEval[solSS, #]) & /@ gCom; eSR = (exteriorValue[rR, #] - analyticEval[solSR, #]) & /@ gCom;
    eTot = (exteriorValue[rS, #] + exteriorValue[rR, #] + lR[#] - analyticEval[chiSol, #]) & /@ gCom;
    ePw = (lPw["f"][#] - analyticEval[solSR, #]) & /@ gCom; eFix = ePw + (gam G0[l, #, r0] & /@ gCom);
    abS = fitPoint9[l, {Re[rS["valR"]] - analyticEval[solSS, aR], Re[rS["derR"]] - analyticEvalD[solSS, aR]}];
    abR = fitPoint9[l, {Re[rR["valR"]] - analyticEval[solSR, aR], Re[rR["derR"]] - analyticEvalD[solSR, aR]}];
    AppendTo[t1Rows, {l, Max[Abs[eSS]]/scale, Abs[abS[[1]]] + Abs[abS[[2]]], Max[Abs[eSR]]/scale, Abs[abR[[1]]] + Abs[abR[[2]]], Max[Abs[eTot]]/scale, Max[Abs[ePw]]/scale, gam}];
    csCheck["Test 1 l=" <> ToString[l] <> ": SS exterior = (1/2) d_z^2 (q^2 Log rho) modes (scaled)", Max[Abs[eSS]]/scale, 10^-10];
    csCheck["Test 1 l=" <> ToString[l] <> ": SS has no point-supported term (|alpha| + |beta|, 1/rho and d_r0(1/rho) modes)", Abs[abS[[1]]] + Abs[abS[[2]]], 10^-12];
    csCheck["Test 1 l=" <> ToString[l] <> ": SR product route exterior (delta carried by the continuation, scaled)", Max[Abs[eSR]]/scale, 10^-10];
    csCheck["Test 1 l=" <> ToString[l] <> ": SR product route has no missing point term (|alpha| + |beta|)", Abs[abR[[1]]] + Abs[abR[[2]]], 10^-12];
    csCheck["Test 1 l=" <> ToString[l] <> ": total exterior = (1/2) d_z^2 chi (scaled)", Max[Abs[eTot]]/scale, 10^-10];
    csCheck["Test 1 l=" <> ToString[l] <> ": SR delta coefficient = -4 Pi q PhiR(x0) Y_l(0)/r0^2", Abs[gam + 4 Pi q (c0 + c1 r0) Y0[l]/r0^2], 10^-30];
    csCheck["Test 1 l=" <> ToString[l] <> ": pointwise real-line SR + delta response = exact (scaled)", Max[Abs[eFix]]/scale, 10^-10];
    Print["    l=", l, ": pointwise SR without the delta misses ", sci[Max[Abs[ePw]]/scale], " (scaled)  [", tick[], " s]"]], {l, lOut9}];

(* ================= TEST 2: Stage 6 puncture ================= *)
Print["\n--- Test 2: Stage 6 puncture q/rho + c q rho^n ---"];
ddR2[n_, c_, L_] := ddR2[n, c, L] = With[{p = aR2[n, c, L]}, Table[D[p, {zVar, k}], {k, 0, 3}]];
ddRes[n_, c_, lp_] := ddResM[n, c, lp, csPar["c2"]];
ddResM[n_, c_, lp_, c2_] := ddResM[n, c, lp, c2] = With[{p = resR[n, c, lp]}, Table[D[p, {zVar, k}], {k, 0, 3}]];
aCdd[n_, c_, L_] := aCdd[n, c, L] = With[{e = aC[n, c, L]}, Table[D[e, {zVar, k}], {k, 0, 2}]];
resCdd[n_, c_, lp_] := resCdd[n, c, lp] = With[{e = resC[n, c, lp]}, Table[D[e, {zVar, k}], {k, 0, 2}]];
(* real mode values (nd derivative orders) and continued (boundary value / complex z) mode values *)
nodeAV[n_, c_, s_, br_, LA_, MV_, nd_: 3] := {Table[pwVal[#, s, br] & /@ Take[ddR2[n, c, L], nd], {L, 0, LA}], Table[pwVal[#, s, br] & /@ Take[ddRes[n, c, lp], nd], {lp, 0, MV}]};
nodeAVc[n_, c_, z_, LA_, MV_] := {Table[evR[#, z] & /@ aCdd[n, c, L], {L, 0, LA}], Table[evR[#, z] & /@ resCdd[n, c, lp], {lp, 0, MV}]};
(* F_j (j <= jMax) for SR = 2 sum_{l' <= M} K^P_{j l'} phi_l' and RR = sum gt phi phi, for every M in Ms *)
nodeF[{A_, V_}, Ms_, rr_: True] := Module[{Mx = Max[Ms], srC, rrSh, tr, z3 = 0 V[[1]]},
   srC = Table[2 prodD[Sum[gt[j, lp, L] A[[L + 1]], {L, Abs[j - lp], j + lp}], V[[lp + 1]]], {j, 0, jMax}, {lp, 0, Mx}];
   tr[j_, l1_, l2_] := With[{g = gt[j, l1, l2]}, If[g == 0, z3, g prodD[V[[l1 + 1]], V[[l2 + 1]]]]];
   rrSh = If[rr, Table[Sum[tr[j, m, l2], {l2, Max[0, m - j], m}] + Sum[tr[j, l1, m], {l1, Max[0, m - j], m - 1}], {j, 0, jMax}, {m, 0, Mx}], None];
   Association@Table[M -> {Total[srC[[All, ;; M + 1]], {2}], If[rr, Total[rrSh[[All, ;; M + 1]], {2}], None]}, {M, Ms}]];
(* Part 1 F_j = 2 sum_{l' <= k} K^P_{j l'} T_l' *)
p1F[A_, tp_, k_, s_] := Table[2 Sum[prodD[Sum[gt[j, lp, L] A[[L + 1]], {L, Abs[j - lp], j + lp}], Take[N[tp[[lp + 1]] /. zVar -> s, 50], Length[A[[1]]]]], {lp, 0, k}], {j, 0, jMax}];
pt[l_, Fl_, s_] := (1/2) dzzArr[Fl, s][[l + 1, 1]];
(* Taylor data of Phi_res at x0 from mode sums: {Phi, d_z, d_z^2, d_x^2, d_z^3, d_z d_x^2} *)
taylor9[n_, c_, M_] := Module[{d},
   d[l_, k_] := evR[D[resR[n, c, l][[1]], {zVar, k}], r0];
   {Sum[Y0[l] d[l, 0], {l, 0, M}], Sum[Y0[l] d[l, 1], {l, 0, M}], Sum[Y0[l] d[l, 2], {l, 0, M}],
    Sum[Y0[l] (d[l, 1]/r0 - l (l + 1)/2 d[l, 0]/r0^2), {l, 0, M}], Sum[Y0[l] d[l, 3], {l, 0, M}],
    Sum[Y0[l] ((d[l, 2]/r0 - d[l, 1]/r0^2) - l (l + 1)/2 (d[l, 1]/r0^2 - 2 d[l, 0]/r0^3)), {l, 0, M}]}];
tModes9[k_, {A0_, A1_, A2_, B2_, A3_, B3_}] := Module[{mu, cl, poly},
   poly = A0 + If[k >= 1, A1 (zVar mu - r0), 0] + If[k >= 2, A2/2 (zVar mu - r0)^2 + B2/2 zVar^2 (1 - mu^2), 0] +
     If[k >= 3, A3/6 (zVar mu - r0)^3 + B3/2 zVar^2 (1 - mu^2) (zVar mu - r0), 0];
   cl = CoefficientList[Expand[poly], mu];
   Table[Expand[Sum[cl[[m + 1]] muPowerMode[lp, m], {m, 0, Length[cl] - 1}]], {lp, 0, k}]];
tPoly9[k_, tc_] := Table[With[{t = tModes9[k, tc][[lp + 1]]}, Table[D[t, {zVar, m}], {m, 0, 3}]], {lp, 0, k}];
(* exact SR + RR Product source: [2 q PhiR/rho] + [PhiR^2] - 2 c q^2 [rho^(n-1)] - c^2 q^2 [rho^(2n)] *)
exPairF[n_, c_] := Function[j, If[j < 0, {0, 0}, pieceAnalyticSource["Product", "SR", j] + pieceAnalyticSource["Product", "RR", j] - 2 c q^2 rhoN[n - 1, j] - c^2 q^2 rhoN[2 n, j]]];
(* product route: shadow term on the axis, Re(continued source)(x + i0) - pointwise, propagated to aR *)
bFrom[l_, M_, nodes_, FR_, FC_] := Module[{s = nodes[[All, 1]], w = nodes[[All, 2]], fb, g0v, g1v},
   fb = Table[Re[pt[l, FC[[i]][M][[1]], s[[i]]]] - pt[l, FR[[i]][M][[1]], s[[i]]], {i, Length[s]}];
   g0v = (s^2/(2 l + 1)) (aR^l s^(-l - 1) - s^l aR^(-l - 1)); g1v = (s^2/(2 l + 1)) (l aR^(l - 1) s^(-l - 1) + (l + 1) s^l aR^(-l - 2));
   {Total[w g0v fb], Total[w g1v fb]}];
resp[l_, d_] := Join[ConstantArray[0, Length[gIn]], homOnGrid[l, aR, d, gOut]];

rows9 = {}; gRows = {}; tRows9 = {}; jRows = {};
Do[Module[{c = r0^-(n + 1), nodesC = cNodes[12], LA = jMax + Mm, FN, AVin, AVout, FIn, FOut, tubeN, FtR, FtC, tcA, ep},
   FN = Table[nodeF[nodeAV[n, c, s, Automatic, LA, Mm], lpList9], {s, nodesC[[All, 1]]}];
   AVin = nodeAV[n, c, r0, "in", LA, Mm, 4]; AVout = nodeAV[n, c, r0, "out", LA, Mm, 4];
   FIn = nodeF[AVin, lpList9]; FOut = nodeF[AVout, lpList9];
   Print["    n=", n, ": ", Length[nodesC], " nodes cached [", tick[], " s]"];
   tubeN = Join[gl[aL, r0, 40], gl[r0, aR, 40]];
   FtR = Table[nodeF[nodeAV[n, c, s, Automatic, jMax + Mprod, Mprod], lpProd, False], {s, tubeN[[All, 1]]}];
   FtC = Table[nodeF[nodeAVc[n, c, s, jMax + Mprod, Mprod], lpProd, False], {s, tubeN[[All, 1]]}];
   Print["    n=", n, ": tube boundary values cached [", tick[], " s]"];
   tcA = Association@Table[M -> taylor9[n, c, M], {M, lpList9}];
   Do[AppendTo[tRows9, {n, M, Sequence @@ Abs[tcA[M] - N[{c0 + c1 r0, c1, 0, 0, 0, 0}, 40]]}], {M, lpList9}];
   ep = exPairF[n, c];
   Do[Module[{chiSol = Sum[pieceAnalyticSolution["Dzz", p, l], {p, {"SS", "SR", "RR"}}], scale, exE, exPt, gEx, g0},
      scale = Max[Abs[analyticEval[chiSol, #] & /@ gCom]];
      exE = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, ep[j][[b]]]], l], {b, 2}];
      exPt = pwVal[exE, #] & /@ nodesC[[All, 1]];
      gEx = (1/2) N[dzzDeltaCoeff[ep, l], 40];
      g0 = G0[l, #, r0] & /@ gCom;
      csCheck["n=" <> ToString[n] <> " l=" <> ToString[l] <> ": exact SR+RR delta coefficient = -4 Pi q PhiR(x0) Y_l(0)/r0^2", Abs[gEx + 4 Pi q (c0 + c1 r0) Y0[l]/r0^2], 10^-30];
      Do[Module[{trunc, gSR, gRR, eC, tc = tcA[M], kRes = {}, jIn, jOut},
         trunc = Table[pt[l, FN[[i]][M][[1]], nodesC[[i, 1]]] + pt[l, FN[[i]][M][[2]], nodesC[[i, 1]]], {i, Length[nodesC]}];
         gSR = (1/2) dzzDeltaNum[l, FIn[M][[1, # + 1, 2]] &, FOut[M][[1, # + 1, 2]] &];
         (* RR: delta coefficient from the one-sided derivatives of the truncated RR at r0 *)
         gRR = (1/2) dzzDeltaNum[l, FIn[M][[2, # + 1, 2]] &, FOut[M][[2, # + 1, 2]] &];
         (* (c): distributional real-line solve of the truncated SR + RR minus exact, zero data at r = 1 *)
         eC = cResponse[l, nodesC, trunc - exPt, gCom] + (gSR + gRR - gEx) g0;
         csCheck["n=" <> ToString[n] <> " M=" <> ToString[M] <> " l=" <> ToString[l] <> ": truncated SR delta = -4 Pi q Y_l(0) PhiRes_M(x0)/r0^2 (mode sum)",
          Abs[gSR + 4 Pi q Y0[l] tc[[1]]/r0^2]/Abs[gEx], 10^-25];
         If[M <= Mprod, Module[{bh = resp[l, bFrom[l, M, tubeN, FtR, FtC]]},
            AppendTo[rows9, {n, "product", M, l, Max[Abs[eC + bh]]/scale, Max[Abs[bh]]/scale, Max[Abs[eC]]/scale}]]];
         (* H0: the Cauchy continuation of the pointwise real-line SR has no delta: the solve misses gSR delta *)
         AppendTo[rows9, {n, "H0", M, l, Max[Abs[eC - gSR g0]]/scale, Max[Abs[gSR g0]]/scale, Max[Abs[eC]]/scale}];
         (* Part 2 one-sided values and derivatives at r0 (H0: Part 1 = 0) *)
         jIn = (1/2) dzzArr[FIn[M][[1]], r0][[l + 1, 1 ;; 2]]; jOut = (1/2) dzzArr[FOut[M][[1]], r0][[l + 1, 1 ;; 2]];
         AppendTo[jRows, {n, M, l, "H0", Abs[jOut[[1]] - jIn[[1]]], Abs[jOut[[2]] - jIn[[2]]]}];
         (* H1(k): Part 1 (closed form, continued) keeps its delta gP1; Part 2 (pointwise + Cauchy) misses gSR - gP1 *)
         Do[Module[{tp = tPoly9[k, tc], pIn, pOut, gP1, d2, pj},
            pIn = p1F[AVin[[1]], tp, k, r0]; pOut = p1F[AVout[[1]], tp, k, r0];
            gP1 = (1/2) dzzDeltaNum[l, pIn[[# + 1, 2]] &, pOut[[# + 1, 2]] &];
            d2 = gSR - gP1;
            AppendTo[rows9, {n, "H1 k=" <> ToString[k], M, l, Max[Abs[eC - d2 g0]]/scale, Max[Abs[d2 g0]]/scale, Max[Abs[eC]]/scale}];
            AppendTo[kRes, Abs[d2]];
            pj = (jOut - jIn) - ((1/2) dzzArr[pOut, r0] - (1/2) dzzArr[pIn, r0])[[l + 1, 1 ;; 2]];
            AppendTo[jRows, {n, M, l, "H1 k=" <> ToString[k], Abs[pj[[1]]], Abs[pj[[2]]]}]], {k, kList9}];
         AppendTo[gRows, {n, M, l, gEx, gSR, gRR, Max[kRes]}]], {M, lpList9}]], {l, lOut9}];
   Print["    n=", n, " done [", tick[], " s]"]], {n, nS9}];
Do[csCheck["n=" <> ToString[r[[1]]] <> " M=" <> ToString[r[[2]]] <> " l=" <> ToString[r[[3]]] <> ": H1 Part 2 delta coefficient |gSR - gP1| (all k) relative to |gSR|", r[[7]]/Abs[r[[5]]], 10^-25], {r, gRows}];
Do[csCheck["n=" <> ToString[r[[1]]] <> " M=" <> ToString[r[[2]]] <> " l=" <> ToString[r[[3]]] <> ": RR has no delta term (relative to |gSR|)", Abs[r[[6]]]/Abs[r[[5]]], 10^-25], {r, Select[gRows, #[[2]] == Mm &]}];

(* ================= validation by contour solves (n = 3, l = 2) ================= *)
Print["\n--- validation of the Test-2 decomposition by contour solves (n = 3, l = 2) ---"];
valRows = {};
Module[{n = 3, c = r0^-4, M = 8, l = 2, LA, tc, tCh, aN, yA, srAt, sCh, sG, sL, AVin, AVout, FIn, FOut, gSR, lineAR, ref, chk},
  LA = jMax + M; tc = taylor9[n, c, M];
  srAt[s_] := pt[l, nodeF[nodeAV[n, c, s, Automatic, LA, M], {M}, False][M][[1]], s];
  AVin = nodeAV[n, c, r0, "in", LA, M]; AVout = nodeAV[n, c, r0, "out", LA, M];
  FIn = nodeF[AVin, {M}, False][M][[1]]; FOut = nodeF[AVout, {M}, False][M][[1]];
  gSR = (1/2) dzzDeltaNum[l, FIn[[# + 1, 2]] &, FOut[[# + 1, 2]] &];
  lineAR = gSR {G0[l, aR, r0], dG0[l, aR, r0]}; ref = Max[Abs[lineAR]];               (* delta part of the distributional real-line SR at aR *)
  tCh = {chebX[cI[[1]], r0, NC], chebX[r0, cI[[2]], NC]}; aN = aCheb[NA]; yA = N[yC[#], 50] & /@ aN;
  sCh = Map[srAt, tCh, {2}]; sG = Map[srAt, gaussP[[All, All, 1]], {2}]; sL = srAt /@ lineNodes[[All, 1]];
  chk[var_, d_, pred_, tol_] := (AppendTo[valRows, {var, M, d, pred}];
    csCheck["validation " <> var <> " (n=3, l=2): contour - distributional line at aR = predicted (relative to |gSR G0|)", Max[Abs[d - pred]]/ref, tol];
    Print["    ", var, ": contour - line = ", sci /@ d, ",  predicted ", sci /@ pred, "  [", tick[], " s]"]);
  Do[Module[{k = If[var === "H0", -1, ToExpression[StringTake[var, -1]]], tp, p1At, p1Ch, p1G, src, gP1},
     If[k < 0, src = cauchyAt[sG, derivCoefs /@ sCh, #] & /@ yA; gP1 = 0,
      tp = tPoly9[k, tc];
      p1At[s_] := pt[l, p1F[nodeAV[n, c, s, Automatic, LA, 0][[1]], tp, k, s], s];
      p1Ch = Map[p1At, tCh, {2}]; p1G = Map[p1At, gaussP[[All, All, 1]], {2}];
      Module[{p1c = (1/2) dzzExpr[Function[j, If[j < 0, 0, 2 Sum[KPc[n, c, j, lp] tp[[lp + 1, 1]], {lp, 0, k}]]], l]},
       src = (evR[p1c, #] & /@ yA) + (cauchyAt[sG - p1G, derivCoefs /@ (sCh - p1Ch), #] & /@ yA)];
      gP1 = (1/2) dzzDeltaNum[l, p1F[AVin[[1]], tp, k, r0][[# + 1, 2]] &, p1F[AVout[[1]], tp, k, r0][[# + 1, 2]] &]];
     chk[var, deltaAt[l, src, sL] - lineAR, -(gSR - gP1) {G0[l, aR, r0], dG0[l, aR, r0]}, 10^-8]], {var, {"H0", "H1 k=0", "H1 k=2"}}];
  (* product route, l'max = 4: continued K~P phi~res on the contour; prediction = shadow term on the axis *)
  Module[{M4 = 4, LA4 = jMax + 4, src, sL4, F4in, F4out, g4, tubeN = Join[gl[aL, r0, 40], gl[r0, aR, 40]], bb, d},
   src = Table[pt[l, nodeF[nodeAVc[n, c, z, LA4, M4], {M4}, False][M4][[1]], z], {z, yA}];
   sL4 = Table[pt[l, nodeF[nodeAV[n, c, s, Automatic, LA4, M4], {M4}, False][M4][[1]], s], {s, lineNodes[[All, 1]]}];
   F4in = nodeF[nodeAV[n, c, r0, "in", LA4, M4], {M4}, False][M4][[1]]; F4out = nodeF[nodeAV[n, c, r0, "out", LA4, M4], {M4}, False][M4][[1]];
   g4 = (1/2) dzzDeltaNum[l, F4in[[# + 1, 2]] &, F4out[[# + 1, 2]] &];
   bb = bFrom[l, M4, tubeN, Table[nodeF[nodeAV[n, c, s, Automatic, LA4, M4], {M4}, False], {s, tubeN[[All, 1]]}],
     Table[nodeF[nodeAVc[n, c, s, LA4, M4], {M4}, False], {s, tubeN[[All, 1]]}]];
   d = deltaAt[l, src, sL4] - g4 {G0[l, aR, r0], dG0[l, aR, r0]};
   AppendTo[valRows, {"product", M4, d, bb}];
   csCheck["validation product route (n=3, l=2, l'max=4): contour - distributional line at aR = axis shadow term (relative)", Max[Abs[d - bb]]/(Max[Abs[bb]] + 10^-10 ref), 10^-6];
   Print["    product: contour - line = ", sci /@ d, ",  predicted ", sci /@ bb, "  [", tick[], " s]"]]];

(* ================= tables and figure ================= *)
mx9[n_, v_, M_, k_] := Max[Cases[rows9, {n, v, M, __}][[All, k]]];
vars9 = Join[{"product", "H0"}, ("H1 k=" <> ToString[#]) & /@ kList9];
md1 = StringRiffle[Join[{"| l | SS (scaled) | SS: abs(alpha)+abs(beta) | SR product route | SR: abs(alpha)+abs(beta) | total | pointwise SR, no delta | delta coefficient |", "|---|---|---|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> StringRiffle[sci /@ #[[2 ;; 7]], " | "] <> " | " <> ToString[NumberForm[N[#[[8]], 6], 6]] <> " |") & /@ t1Rows], "\n"];
md2 = StringRiffle[Join[{"| n | route | " <> StringRiffle[ToString /@ lpList9, " | "] <> " |", "|---|---|" <> StringRepeat["---|", Length[lpList9]]},
    Flatten@Table[Join[{"| " <> ToString[n] <> " | (c) | " <> StringRiffle[sci[mx9[n, "H0", #, 7]] & /@ lpList9, " | "] <> " |"},
       Table["| " <> ToString[n] <> " | " <> v <> " | " <> StringRiffle[If[v === "product" && # > Mprod, "-", sci[mx9[n, v, #, 5]]] & /@ lpList9, " | "] <> " |", {v, vars9}]], {n, nS9}]], "\n"];
md3 = StringRiffle[Join[{"| n | l | exact SR+RR | truncated SR (l'max = " <> ToString[Mm] <> ") | truncated RR | Part 2 of H1, max over k |", "|---|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> ToString[#[[3]]] <> " | " <> ToString[NumberForm[N[#[[4]], 8], 8]] <> " | " <> ToString[NumberForm[N[#[[5]], 8], 8]] <> " | " <> sci[#[[6]]] <> " | " <> sci[#[[7]]] <> " |") & /@ Select[gRows, #[[2]] == Mm &]], "\n"];
md4 = StringRiffle[Join[{"| n | route | max_l jump of Part 2 at r0 | max_l jump of its r-derivative |", "|---|---|---|---|"},
    Flatten@Table["| " <> ToString[n] <> " | " <> v <> " | " <> sci[Max[Cases[jRows, {n, Mm, _, v, a_, _} :> a]]] <> " | " <> sci[Max[Cases[jRows, {n, Mm, _, v, _, b_} :> b]]] <> " |", {n, nS9}, {v, Rest[vars9]}]], "\n"];
md5 = StringRiffle[Join[{"| route | l'max | contour - line (value, derivative) at aR | predicted |", "|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3]], ", "] <> " | " <> StringRiffle[sci /@ #[[4]], ", "] <> " |") & /@ valRows], "\n"];
md6 = StringRiffle[Join[{"| n | l'max | Phi(x0) | d_z | d_z^2 | d_x^2 | d_z^3 | d_z d_x^2 |", "|---|---|---|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;;]], " | "] <> " |") & /@ tRows9], "\n"];
Print["\nTest 1 (exact puncture):\n", md1, "\n\nTest 2 (max over l of scaled exterior error):\n", md2, "\n\ndelta coefficients at l'max = ", Mm, ":\n", md3,
  "\n\nPart 2 smoothness at r0 (l'max = ", Mm, "):\n", md4, "\n\nvalidation:\n", md5, "\n\nTaylor coefficient errors (mode sums):\n", md6];
Export[FileNameJoin[{resDir, "stage9_summary.md"}],
  "# Stage 9: Q = (1/2) d_z^2 (A B)\n\n## Test 1: exact puncture, product route (scaled exterior errors)\n\n" <> md1 <>
   "\n\n## Test 2: Stage 6 puncture, R = 1/2 (max over l of scaled exterior error; (c) = l' truncation)\n\n" <> md2 <>
   "\n\n## delta(r - r0) coefficients of the source modes\n\n" <> md3 <> "\n\n## Part 2 smoothness at r0\n\n" <> md4 <>
   "\n\n## validation by contour solves (n = 3, l = 2)\n\n" <> md5 <> "\n\n## Taylor coefficient errors (mode sums)\n\n" <> md6 <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage9_test1.csv"}], Prepend[{#[[1]], Sequence @@ (sci /@ #[[2 ;; 7]]), N[#[[8]], 12]} & /@ t1Rows,
   {"l", "SS_err", "SS_alpha_beta", "SR_err", "SR_alpha_beta", "total_err", "pointwise_no_delta_err", "delta_coeff"}], "CSV"];
Export[FileNameJoin[{resDir, "stage9_test2.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]]} & /@ rows9, {"n", "route", "lpmax", "l", "total", "route_specific", "c"}], "CSV"];
Export[FileNameJoin[{resDir, "stage9_delta.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], N[#[[4]], 12], N[#[[5]], 12], sci[#[[6]]], sci[#[[7]]]} & /@ gRows, {"n", "lpmax", "l", "gamma_exact", "gamma_SR_trunc", "gamma_RR_trunc", "gamma_part2_maxk"}], "CSV"];
Export[FileNameJoin[{resDir, "stage9_part2_jumps.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], sci[#[[5]]], sci[#[[6]]]} & /@ jRows, {"n", "lpmax", "l", "route", "jump", "jump_derivative"}], "CSV"];
Export[FileNameJoin[{resDir, "stage9_taylor.csv"}], Prepend[{#[[1]], #[[2]], Sequence @@ (sci /@ #[[3 ;;]])} & /@ tRows9, {"n", "lpmax", "err_Phi", "err_dz", "err_dzz", "err_dxx", "err_dzzz", "err_dz_dxx"}], "CSV"];
$fs = {FontFamily -> "Times", FontSize -> 9};
vSt = <|"product" -> Directive[Red, AbsoluteThickness[1.4]], "H0" -> Directive[Black, AbsoluteThickness[1.2]], "H1 k=0" -> Directive[RGBColor[0.55, 0.3, 0.7], AbsoluteThickness[2.6]],
   "H1 k=1" -> Directive[RGBColor[0.85, 0.45, 0.1], AbsoluteThickness[1.9]], "H1 k=2" -> Directive[RGBColor[0.1, 0.35, 0.75], AbsoluteThickness[1.3]], "H1 k=3" -> Directive[RGBColor[0.15, 0.6, 0.25], AbsoluteThickness[0.8]]|>;
panels = Table[ListLogPlot[Join[Table[Table[{M, mx9[n, v, M, 5]}, {M, If[v === "product", lpProd, lpList9]}], {v, vars9}], {Table[{M, mx9[n, "H0", M, 7]}, {M, lpList9}]}],
    Joined -> True, PlotMarkers -> Join[Table[{"\[FilledCircle]", 5}, {Length[vars9]}], {None}],
    PlotStyle -> Join[Values[vSt], {Directive[Gray, Dashed, AbsoluteThickness[1]]}], PlotRange -> All,
    Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]'\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 212,
    Epilog -> {Inset[Style["n = " <> ToString[n], $fs], Scaled[{0.05, 0.05}], {Left, Bottom}]}], {n, nS9}];
leg = LineLegend[Join[Values[vSt], {Directive[Gray, Dashed]}], {"product route", "H0 (no split)", "H1, k = 0", "H1, k = 1", "H1, k = 2", "H1, k = 3", "(c) truncation"}, LegendLayout -> {"Row", 4}, LabelStyle -> $fs];
Export[FileNameJoin[{figDir, "fig_dzz.pdf"}], Legended[Grid[{panels}, Spacings -> 1], Placed[leg, Below]]];
Export[FileNameJoin[{figDir, "fig_dzz.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], N[#[[5]], 6], N[#[[6]], 6], N[#[[7]], 6]} & /@ rows9, {"n", "route", "lpmax", "l", "total", "route_specific", "c"}], "CSV"];
Print["\nStage 9 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage9_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
