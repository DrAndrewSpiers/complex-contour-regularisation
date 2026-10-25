#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 14: cost of the SS angular integrals, real axis against contour (Dzz case).
   S(r, mu) = (1/2) q^2 d_z^2 rho^-2 = (1/2) q^2 (-2/rho^4 + 8 (r mu - r0)^2/rho^6);  SS_l(r) = Sqrt[Pi (2l+1)] Int S P_l dmu,  l = 0..Lmax.
   Every evaluation of the 3D integrand S is counted (sCount); P_l evaluations are not.
   (A) real axis, adaptive NIntegrate in mu per (l, r), PrecisionGoal 8 and 12, r = r0 +- k h, h = R/200, k = 1..200.
   (B) contour, Stage 12 (D) machinery: rho path z + r0 -> -+i/2 -> r0 - z (in) or z - r0 (out), nq Gauss nodes per segment, each l separately.
   (C) contour, shared nodes: S once per rho node, all l by the Legendre recurrence.
   (D) real axis, shared nodes: Gauss in u = ln rho from ln|Dr| to ln(r + r0) (and plain rho for comparison).
   (E) radial nodes: Gauss panels in a on the contour needed for the exterior SS field (Green's-function quadrature).
   Errors: pointwise relative to the closed form (40-digit evaluation), max over l <= Lmax and points.
   Run:  wolframscript -file stage14_cost.wls                                                                          *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]]; Needs["NumericalDifferentialEquationAnalysis`"];
resDir = FileNameJoin[{dir, "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>];
q = 1; r0 = 5; R = 1/2; aL = r0 - R; aR = r0 + R;
csResetChecks[]; Off[General::munfl]; t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
Lt = 30; LmaxList = {10, 20, 30}; tols = {10^-8, 10^-12};
nqList = {8, 12, 16, 20, 24, 32, 40, 48, 56, 64, 72, 80, 96, 112, 128, 160};
$fs = {FontFamily -> "Times", FontSize -> 9};
Print["=== Stage 14: SS angular-integral cost, real axis vs contour; Lmax = ", LmaxList, " ==="];

(* ---------------- counted 3D integrand ---------------- *)
sCount = 0;
SmuN[r_, mu_?NumericQ] := (sCount++; With[{rho2 = r^2 + r0^2 - 2 r r0 mu}, (1/2) q^2 (-2/rho2^2 + 8 (r mu - r0)^2/rho2^3)]);
Srho[z_, rho_List] := (sCount += Length[rho]; With[{d = (z^2 - r0^2 - rho^2)/(2 r0)}, (1/2) q^2 (-2/rho^4 + 8 d^2/rho^6)]);
muOf[z_, rho_] := (z^2 + r0^2 - rho^2)/(2 z r0);
legAll[L_, mu_] := Module[{P = {ConstantArray[1., Length[mu]], mu}}, Do[AppendTo[P, ((2 l + 1) mu P[[l + 1]] - l P[[l]])/(l + 1)], {l, 1, L - 1}]; P];
gq[n_] := gq[n] = Transpose[N[GaussianQuadratureWeights[n, 0, 1, 30]]];
nrm = N[Table[Sqrt[Pi (2 l + 1)], {l, 0, Lt}]];

(* ---------------- closed-form references ---------------- *)
sRealE[l_] := ssSourceModeExpr["Dzz", l, zVar, "real"] /. Log[x_^2] :> 2 Log[Abs[x]];
refRe[l_] := refRe[l] = Function @@ {{zVar}, N[sRealE[l], 60]};
refOut[l_] := refOut[l] = Function @@ {{zVar}, N[ssSourceModeExpr["Dzz", l, zVar, "continued"], 60]};        (* continuation from a_R *)
refIn[l_] := refIn[l] = Function @@ {{zVar}, N[sRealE[l] /. Abs[x_] :> -x, 60]};                            (* continuation from a_L *)
ev[f_, z_] := N[f[SetPrecision[z, 40]]];
Do[csCheck["l=" <> ToString[l] <> ": in/out references agree with the real modes on their sides (r = 4.7, 5.3)",
   Abs[ev[refIn[l], 47/10] - ev[refRe[l], 47/10]]/Abs[ev[refRe[l], 47/10]] + Abs[ev[refOut[l], 53/10] - ev[refRe[l], 53/10]]/Abs[ev[refRe[l], 53/10]], 10^-14], {l, {0, 10, 30}}];
relErr[v_, ref_] := Abs[v - ref]/Abs[ref];

(* ---------------- contour nodes (Stage 12) ---------------- *)
yA[a_] := a + I (1 - ((a - r0)/R)^2)^6;
aNodes[P_, m_] := Module[{br = Subdivide[aL, aR, P], a}, Flatten[Table[({N[yA[#[[1]]]], N[D[yA[a], a] /. a -> #[[1]]], #[[2]]} & /@ N[GaussianQuadratureWeights[m, br[[k]], br[[k + 1]], 30]]), {k, P}], 1]];
zAll = aNodes[40, 30][[All, 1]]; zS = zAll[[;; ;; 4]];                                                      (* 300 sample nodes for accuracy *)
refC = Association@Table[br -> Table[ev[If[br === "in", refIn[l], refOut[l]], z], {z, zS}, {l, 0, Lt}], {br, {"in", "out"}}];
Print["    references on the contour [", tick[], " s]"];

(* path nodes and weights for one z and branch (two segments of nq Gauss nodes) *)
pathNodes[z_, br_, nq_] := Module[{via = If[br === "in", -I/2, I/2], e = If[br === "in", r0 - z, z - r0], t, w}, {t, w} = gq[nq];
   {Join[z + r0 + (via - z - r0) t, via + (e - via) t], Join[(via - z - r0) w, (e - via) w]}];
(* (C) shared: all l from one set of S evaluations *)
modesC[z_, br_, nq_, L_] := Module[{rho, w, s, mu}, {rho, w} = pathNodes[z, br, nq]; s = Srho[z, rho]; mu = muOf[z, rho];
   nrm[[;; L + 1]] (legAll[L, mu] . (w s (-rho/(z r0))))];
(* (B) per l: S recomputed for every l, LegendreP built in *)
modeB[l_, z_, br_, nq_] := Module[{rho, w, s, mu}, {rho, w} = pathNodes[z, br, nq]; s = Srho[z, rho]; mu = muOf[z, rho];
   nrm[[l + 1]] Total[w s LegendreP[l, mu] (-rho/(z r0))]];

(* ================= (B), (C): accuracy against nq on the contour ================= *)
errC = Association@Table[nq -> Module[{e = Table[0., {Lt + 1}]},
      Do[Do[Module[{v = modesC[zS[[i]], br, nq, Lt]}, e = MapThread[Max, {e, relErr[v, refC[br][[i]]]}]], {i, Length[zS]}], {br, {"in", "out"}}]; e], {nq, nqList}];
nqReqL[l_, tol_] := SelectFirst[nqList, errC[#][[l + 1]] <= tol &, Infinity];
nqReq[L_, tol_] := SelectFirst[nqList, Max[errC[#][[;; L + 1]]] <= tol &, Infinity];
Print["    contour accuracy scan [", tick[], " s]"];
(* confirm (B) = (C) numerically and the S counts *)
Module[{z = zS[[37]], nq = 80, cB, cC, vB, vC},
  sCount = 0; vB = Table[modeB[l, z, "out", nq], {l, 0, Lt}]; cB = sCount;
  sCount = 0; vC = modesC[z, "out", nq, Lt]; cC = sCount;
  csCheck["(B) per-l and (C) shared nodes give the same modes (nq = 80, l <= 30, relative)", Max[relErr[vB, vC]], 10^-11];
  csCheck["(B) S evaluations per z = 2 nq (Lmax + 1); (C) = 2 nq", Abs[cB - 2 nq (Lt + 1)] + Abs[cC - 2 nq], 0];
  cnt = Table[sCount = 0; modesC[z, "out", nq, L]; sCount, {L, LmaxList}];
  csCheck["(C) S evaluations per z independent of Lmax (nq = 80): " <> ToString[cnt], Max[cnt] - Min[cnt], 0]];

(* ================= (E): radial nodes for the exterior SS field ================= *)
gOut = Subdivide[13/2, 10, 14];
extField[nodes_, nqS_] := Module[{z = nodes[[All, 1]], dzw = nodes[[All, 2]] nodes[[All, 3]], src},
   src = Transpose[modesC[#, "out", nqS, Lt] & /@ z];                                          (* [l][node] *)
   Table[Module[{M1 = Re[Total[z^(1 - l) src[[l + 1]] dzw]], M2 = Re[Total[z^(l + 2) src[[l + 1]] dzw]]}, (#^l M1 - #^(-l - 1) M2)/(2 l + 1) & /@ gOut], {l, 0, Lt}]];
nqE = 128;
refE = extField[aNodes[80, 30], nqE];
eScan = Flatten[Table[Module[{f = extField[aNodes[P, m], nqE]},
      {P, m, P m, Table[Max[Table[Max[Abs[f[[l + 1]] - refE[[l + 1]]]]/Max[Abs[refE[[l + 1]]]], {l, 0, L}]], {L, LmaxList}]}], {P, {1, 2, 4, 6, 8, 12, 16, 24, 32}}, {m, {6, 8, 10, 12, 16, 20, 24}}], 1];
nodesReq[L_, tol_] := Module[{k = Position[LmaxList, L][[1, 1]], ok}, ok = Select[eScan, #[[4, k]] <= tol &]; If[ok === {}, {Infinity, 32, 24}, First[SortBy[ok, #[[3]] &]][[{3, 1, 2}]]]];
csCheck["(E) reference exterior field converged (40x30 vs 80x30 panels, l <= 30, relative)", Max[Table[Max[Abs[(extField[aNodes[40, 30], nqE][[l + 1]] - refE[[l + 1]])]]/Max[Abs[refE[[l + 1]]]], {l, 0, Lt}]], 10^-12];
Print["    radial-node scan [", tick[], " s]"];

(* ================= (D): real axis, shared nodes ================= *)
h = R/200; rPts = Join[r0 - h Range[200], r0 + h Range[200]];
refR = Table[ev[refRe[l], r], {r, rPts}, {l, 0, Lt}];
Print["    real-axis references [", tick[], " s]"];
modesD[r_, nq_, L_, mode_] := Module[{t, w, a = Abs[r - r0], b = r + r0, rho, wr, s, mu},
   {t, w} = gq[nq];
   If[mode === "log", rho = Exp[Log[a] + (Log[b] - Log[a]) t]; wr = (Log[b] - Log[a]) w rho, rho = a + (b - a) t; wr = (b - a) w];
   s = Srho[r, rho]; mu = muOf[r, rho];
   nrm[[;; L + 1]] (legAll[L, mu] . (wr s rho/(r r0)))];
errD = Association@Table[mode -> Association@Table[nq -> Table[relErr[modesD[rPts[[i]], nq, Lt, mode], refR[[i]]], {i, Length[rPts]}], {nq, nqList}], {mode, {"log", "plain"}}];
nqReqD[i_, L_, tol_, mode_] := SelectFirst[nqList, Max[errD[mode][#][[i, ;; L + 1]]] <= tol &, Infinity];
Print["    real-axis shared-node scan [", tick[], " s]"];

(* ================= (A): real axis, adaptive NIntegrate per (l, r) ================= *)
aRes = Association[];
Do[Module[{evals = Table[0, {Length[rPts]}, {Lt + 1}], errs = Table[0., {Length[rPts]}, {Lt + 1}], tm},
   tm = First@AbsoluteTiming[Do[Module[{v},
        sCount = 0; v = Quiet[nrm[[l + 1]] NIntegrate[SmuN[rPts[[i]], mu] LegendreP[l, mu], {mu, -1, 1}, PrecisionGoal -> pg, MaxRecursion -> 50]];
        evals[[i, l + 1]] = sCount; errs[[i, l + 1]] = relErr[v, refR[[i, l + 1]]]], {i, Length[rPts]}, {l, 0, Lt}]];
   aRes[pg] = <|"evals" -> evals, "errs" -> errs, "time" -> tm|>;
   Print["    (A) PrecisionGoal ", pg, ": ", Total[evals, 2], " S evaluations, ", Round[tm, 0.1], " s [", tick[], " s]"]], {pg, {8, 12}}];

(* ================= wall times of the production configurations ================= *)
timeC[L_, tol_] := Module[{nq = Min[nqReq[L, tol], 160], nodes = aNodes @@ Rest[nodesReq[L, tol]]}, First@AbsoluteTiming[modesC[#, "out", nq, L] & /@ nodes[[All, 1]]]];
timeB[L_, tol_] := Module[{nodes = aNodes @@ Rest[nodesReq[L, tol]]}, First@AbsoluteTiming[Do[modeB[l, z, "out", Min[nqReqL[l, tol], 160]], {z, nodes[[All, 1]]}, {l, 0, L}]]];
timeD[L_, tol_] := Module[{nqs = Table[Min[nqReqD[i, L, tol, "log"], 160], {i, Length[rPts]}]}, First@AbsoluteTiming[Do[modesD[rPts[[i]], nqs[[i]], L, "log"], {i, Length[rPts]}]]];

(* ================= summary ================= *)
rows = {};
Do[Module[{k = Position[LmaxList, L][[1, 1]]},
   Do[Module[{pg = If[tol == 10^-8, 8, 12], ev, er},
      ev = Total[aRes[pg]["evals"][[All, ;; L + 1]], 2]; er = Max[aRes[pg]["errs"][[All, ;; L + 1]]];
      AppendTo[rows, {L, sci[N[tol]], "(A) real axis, adaptive per (l, r), PrecisionGoal " <> ToString[pg], Length[rPts], ev, sci[er], N[aRes[pg]["time"] (L + 1)/(Lt + 1), 3]}]], {tol, tols}];
   Do[Module[{nqs = Table[nqReqD[i, L, tol, "log"], {i, Length[rPts]}], ev},
      ev = If[MemberQ[nqs, Infinity], Infinity, Total[nqs]];
      AppendTo[rows, {L, sci[N[tol]], "(D) real axis, shared nodes in ln rho", Length[rPts], ev,
        sci[Max[Table[Max[errD["log"][Min[nqs[[i]], 160]][[i, ;; L + 1]]], {i, Length[rPts]}]]], N[timeD[L, tol], 3]}]], {tol, tols}];
   Do[Module[{ne = nodesReq[L, tol]},
      AppendTo[rows, {L, sci[N[tol]], "(B) contour, per l (2 nq per segment pair)", ne[[1]], ne[[1]] Total[Table[2 nqReqL[l, tol], {l, 0, L}]],
        sci[Max[Table[errC[Min[nqReqL[l, tol], 160]][[l + 1]], {l, 0, L}]]], N[timeB[L, tol], 3]}];
      AppendTo[rows, {L, sci[N[tol]], "(C) contour, shared nodes", ne[[1]], ne[[1]] 2 nqReq[L, tol], sci[Max[errC[Min[nqReq[L, tol], 160]][[;; L + 1]]]], N[timeC[L, tol], 3]}]], {tol, tols}]], {L, LmaxList}];
csCheck["(C) reaches 1e-12 for Lmax = 20 (nq <= 160)", If[nqReq[20, 10^-12] < Infinity, 0, 1], 0];
csCheck["(C) Lmax = 30: machine-precision floor of the shared-node quadrature, min over nq of the max error (reported, < 3e-12)", Min[Table[Max[errC[nq]], {nq, nqList}]], 3 10^-12];
md1 = StringRiffle[Join[{"| Lmax | target | method | radial points | S evaluations | max error | wall time (s) |", "|---|---|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> #[[2]] <> " | " <> #[[3]] <> " | " <> ToString[#[[4]]] <> " | " <> ToString[#[[5]]] <> " | " <> #[[6]] <> " | " <> ToString[#[[7]]] <> " |") & /@ rows], "\n"];
(* per-point cost against |Dr| *)
dSel = {1, 2, 5, 10, 20, 50, 100, 200};
md2 = StringRiffle[Join[{"| abs(Dr)/h | (A) evals per (l, r), PG 8 | PG 12 | (D) ln-rho nq, 1e-8 | 1e-12 | (D) plain rho nq, 1e-8 | 1e-12 |", "|---|---|---|---|---|---|---|"},
    Table[Module[{iL = 201 - k, iR = 200 + k},
      "| " <> ToString[k] <> " | " <> ToString[Round[Mean[Flatten[aRes[8]["evals"][[{iL, iR}]]]]]] <> " | " <> ToString[Round[Mean[Flatten[aRes[12]["evals"][[{iL, iR}]]]]]] <> " | " <>
       StringRiffle[ToString /@ {Max[nqReqD[iL, Lt, 10^-8, "log"], nqReqD[iR, Lt, 10^-8, "log"]], Max[nqReqD[iL, Lt, 10^-12, "log"], nqReqD[iR, Lt, 10^-12, "log"]],
          Max[nqReqD[iL, Lt, 10^-8, "plain"], nqReqD[iR, Lt, 10^-8, "plain"]], Max[nqReqD[iL, Lt, 10^-12, "plain"], nqReqD[iR, Lt, 10^-12, "plain"]]}, " | "] <> " |"], {k, dSel}]], "\n"];
md3 = StringRiffle[Join[{"| Lmax | nq per segment for 1e-8 | 1e-12 | contour nodes for 1e-8 (panels x m) | 1e-12 |", "|---|---|---|---|---|"},
    Table["| " <> ToString[L] <> " | " <> ToString[nqReq[L, 10^-8]] <> " | " <> ToString[nqReq[L, 10^-12]] <> " | " <>
      ToString[nodesReq[L, 10^-8][[1]]] <> " (" <> ToString[nodesReq[L, 10^-8][[2]]] <> "x" <> ToString[nodesReq[L, 10^-8][[3]]] <> ") | " <>
      ToString[nodesReq[L, 10^-12][[1]]] <> " (" <> ToString[nodesReq[L, 10^-12][[2]]] <> "x" <> ToString[nodesReq[L, 10^-12][[3]]] <> ") |", {L, LmaxList}]], "\n"];
Print["\n", md1, "\n\n(A)/(D) per point (Lmax = 30):\n", md2, "\n\n(B)/(C)/(E):\n", md3];
Export[FileNameJoin[{resDir, "stage14_summary.md"}], "# Stage 14: SS angular-integral cost, real axis vs contour\n\n" <> md1 <> "\n\n## per-point cost against |Dr| (Lmax = 30, h = R/200)\n\n" <> md2 <>
   "\n\n## contour: nq per segment and radial nodes\n\n" <> md3 <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage14_summary.csv"}], Prepend[rows /. Infinity -> "inf", {"Lmax", "target", "method", "radial_points", "S_evaluations", "max_error", "wall_time_s"}], "CSV"];
Export[FileNameJoin[{resDir, "stage14_contour_nq.csv"}], Prepend[Flatten[Table[{nq, l, sci[errC[nq][[l + 1]]]}, {nq, nqList}, {l, 0, Lt}], 1], {"nq_per_segment", "l", "max_rel_error"}], "CSV"];
Export[FileNameJoin[{resDir, "stage14_radial_nodes.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], Sequence @@ (sci /@ #[[4]])} & /@ eScan, {"panels", "nodes_per_panel", "total_nodes", "err_L10", "err_L20", "err_L30"}], "CSV"];
Export[FileNameJoin[{resDir, "stage14_real_axis.csv"}], Prepend[Flatten[Table[{N[rPts[[i]] - r0, 8], aRes[8]["evals"][[i, l + 1]], aRes[12]["evals"][[i, l + 1]], sci[aRes[8]["errs"][[i, l + 1]]], sci[aRes[12]["errs"][[i, l + 1]]]}, {i, Length[rPts]}, {l, 0, Lt}], 1],
   {"Dr", "evals_PG8", "evals_PG12", "err_PG8", "err_PG12"}], "CSV"];
(* figure: total S evaluations against max error, one panel per Lmax *)
Module[{panels, neC},
  panels = Table[Module[{k = Position[LmaxList, L][[1, 1]], ne = nodesReq[L, 10^-12][[1]], cC, cB, cD, pA},
     cC = Table[{ne 2 nq, Max[errC[nq][[;; L + 1]]]}, {nq, nqList}];
     cB = Table[{ne (L + 1) 2 nq, Max[errC[nq][[;; L + 1]]]}, {nq, nqList}];
     cD = Table[{Length[rPts] nq, Max[Table[Max[errD["log"][nq][[i, ;; L + 1]]], {i, Length[rPts]}]]}, {nq, nqList}];
     pA = Table[{Total[aRes[pg]["evals"][[All, ;; L + 1]], 2], Max[aRes[pg]["errs"][[All, ;; L + 1]]]}, {pg, {8, 12}}];
     ListLogLogPlot[{cC, cB, cD, pA}, Joined -> {True, True, True, False}, PlotMarkers -> {{"\[FilledCircle]", 4}, {"\[FilledSquare]", 4}, {"\[FilledDiamond]", 5}, {"\[FilledUpTriangle]", 8}},
      PlotStyle -> {RGBColor[0.15, 0.6, 0.25], RGBColor[0.1, 0.35, 0.75], RGBColor[0.85, 0.45, 0.1], Red}, Frame -> True, Axes -> False, PlotRange -> {All, {10^-16, 1}},
      FrameLabel -> {"S evaluations", "max relative error"}, BaseStyle -> $fs, ImageSize -> 250, FrameTicks -> {{Automatic, None}, {Table[{10^k, Superscript[10, k]}, {k, 3, 8}], None}},
      Epilog -> {Inset[Style["Lmax = " <> ToString[L] <> " (contour: " <> ToString[ne] <> " nodes)", $fs], Scaled[{0.97, 0.03}], {Right, Bottom}]}]], {L, LmaxList}];
  Export[FileNameJoin[{figDir, "fig_cost.pdf"}], Legended[Grid[{panels}, Spacings -> 0.5],
     Placed[LineLegend[{RGBColor[0.15, 0.6, 0.25], RGBColor[0.1, 0.35, 0.75], RGBColor[0.85, 0.45, 0.1], Red},
       {"(C) contour, shared nodes", "(B) contour, per \[ScriptL]", "(D) real axis, shared nodes (ln \[Rho])", "(A) real axis, adaptive per (\[ScriptL], r), PG 8 and 12"}, LabelStyle -> $fs, LegendLayout -> {"Row", 2}], Below]]];
  Export[FileNameJoin[{figDir, "fig_cost.csv"}], Prepend[rows /. Infinity -> "inf", {"Lmax", "target", "method", "radial_points", "S_evaluations", "max_error", "wall_time_s"}], "CSV"]];
Print["\nStage 14 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage14_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
