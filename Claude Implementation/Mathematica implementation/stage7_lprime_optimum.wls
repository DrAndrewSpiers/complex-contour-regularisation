#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 7: optimum l'max and shrinking worldtube, Stage 6 puncture PhiP = q/rho + c q rho^n, c = r0^-(n+1).
   Total exterior SR+RR error = (c) l'-truncation (real line; independent of R) + (b) shadow contamination (depends on R, l'max).
   Both are computed from the Lemma-1 decomposition by quadrature (validated against Stage 6 three-region/real-line solves here):
     (c): Green's-function quadrature over [1,10] of Delta = (exact SR_P + RR_P source) - (l'-truncated source), zero data at r = 1;
     (b): response at aR of -2 Sum Im K~P Im phi~res (+ gradient analogue) over [r0-R, r0+R].
   R in {1/2, 1/4} and R = r0/l'max; contour amplitude amp = 2R (only used by the three-region verification runs).
   Run:  wolframscript -file stage7_lprime_optimum.wls                                                               *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
s6 = Import[FileNameJoin[{dir, "stage6_3d_puncture_test.wls"}], "Text"];
s6 = StringTake[s6, {1, StringPosition[s6, "(* ================= main"][[1, 1]] - 1}];
s6 = StringReplace[s6, {"dir = DirectoryName[$InputFileName];" -> "", "Print[\"=== Stage 6" -> "Hold[\"x"}];
ToExpression[s6];                                           (* Stage 6 definitions: rhoN, cont, KP*, res*, ssPsrc, runSetting, decomp, ... *)
figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
nS = {3, 5}; lpList = {4, 8, 12, 16, 20, 30, 40}; lpTop = Max[lpList]; lOut = Range[0, 6]; Rfixed = {1/2, 1/4};
smoke = MemberQ[$ScriptCommandLine, "smoke"];       (* wolframscript -file stage7_lprime_optimum.wls smoke *)
If[smoke, nS = {3}; lpList = {4, 8}; lpTop = 8; lOut = {0, 2}];
gIn = Subdivide[1, 7/2, 12]; gOut = Subdivide[13/2, 10, 14]; gCom = Join[gIn, gOut];        (* common exterior grid, valid for R <= 3/2 *)
Print["=== Stage 7: l'max optimum and shrinking tube; n = ", nS, ", l'max = ", lpList, ", l = ", lOut, " ==="];

(* ---------------- node values (all real-axis quantities) ---------------- *)
gl[a_, b_, m_] := GaussianQuadratureWeights[m, a, b, 50];
cNodes[m_, grid_: gCom] := Module[{bp = Union[Join[{1, 4, 9/2, 11/2, 6, r0, 10}, grid]]}, Join @@ Table[gl[bp[[k]], bp[[k + 1]], m], {k, Length[bp] - 1}]];
modeVals[Q_, n_, c_, s_] := Module[{Ls = Range[0, 6 + Max[lpTop, 12]], lps = Range[0, Max[lpTop, 12]]},
   <|"a" -> (pwVal[aR2[n, c, #], s] & /@ Ls), "da" -> (pwVal[daR2[n, c, #], s] & /@ Ls),
     "v" -> (pwVal[resR[n, c, #], s] & /@ lps), "dv" -> (pwVal[dresR[n, c, #], s] & /@ lps)|>];
daR2[n_, c_, L_] := daR2[n, c, L] = D[#, zVar] & /@ aR2[n, c, L];
angF[l_, l1_, l2_] := (l1 (l1 + 1) + l2 (l2 + 1) - l (l + 1))/2;
(* truncated SR_P + RR_P source at one node from cached mode values *)
truncSrc[Q_, l_, M_, mv_, s_] := Module[{a = mv["a"], da = mv["da"], v = mv["v"], dv = mv["dv"], sr, rr},
   sr = 2 Sum[Module[{Ls = Range[Abs[l - lp], l + lp]},
       If[Q === "Product", Sum[gt[l, lp, L] a[[L + 1]], {L, Ls}] v[[lp + 1]],
        Sum[gt[l, lp, L] da[[L + 1]], {L, Ls}] dv[[lp + 1]] + Sum[gt[l, lp, L] angF[l, lp, L] a[[L + 1]], {L, Ls}] v[[lp + 1]]/s^2]], {lp, 0, M}];
   rr = Sum[If[gt[l, l1, l2] == 0, 0, gt[l, l1, l2] If[Q === "Product", v[[l1 + 1]] v[[l2 + 1]],
        dv[[l1 + 1]] dv[[l2 + 1]] + angF[l, l1, l2] v[[l1 + 1]] v[[l2 + 1]]/s^2]], {l1, 0, M}, {l2, Max[0, Abs[l - l1]], Min[M, l + l1]}];
   sr + rr];
(* (c)(r) on a grid: error e with L e = (truncated - exact source), zero data at r = 1,
   e(r) = Int_1^r s^2/(2l+1) (r^l s^(-l-1) - s^l r^(-l-1)) delta(s) ds, by cumulative Gauss quadrature (grid points are breakpoints) *)
cResponse[l_, nodes_, delta_, grid_] := Table[Module[{sel = Select[Transpose[{nodes[[All, 1]], nodes[[All, 2]], delta}], #[[1]] < r &]},
     Total[(#[[2]] #[[1]]^2/(2 l + 1) (r^l #[[1]]^(-l - 1) - #[[1]]^l r^(-l - 1)) #[[3]]) & /@ sel]], {r, grid}];

(* (b): per-l' responses at aR for a tube of radius R (fast version of Stage 6 decomp, cached node values) *)
bPerLp[Q_, n_, c_, R_, l_, lpMax_, np_] := Module[{aLv = r0 - R, aRv = r0 + R, nd, s, w, g0, g1, aCv, daCv, rc, drc, Ls},
   nd = Join[gl[aLv, r0, np], gl[r0, aRv, np]]; s = nd[[All, 1]]; w = nd[[All, 2]];
   Ls = Range[0, l + lpMax];
   aCv = Table[evR[aC[n, c, L], #] & /@ s, {L, Ls}]; daCv = If[Q === "Product", None, Table[evR[daC[n, c, L], #] & /@ s, {L, Ls}]];
   g0 = (s^2/(2 l + 1)) (aRv^l s^(-l - 1) - s^l aRv^(-l - 1));
   g1 = (s^2/(2 l + 1)) (l aRv^(l - 1) s^(-l - 1) + (l + 1) s^l aRv^(-l - 2));
   Table[Module[{rcv = evR[resC[n, c, lp], #] & /@ s, kc, fb},
      kc = Sum[gt[l, lp, L] aCv[[L + 1]], {L, Range[Abs[l - lp], l + lp]}];
      fb = If[Q === "Product", -2 Im[kc] Im[rcv],
        Module[{drcv = evR[D[resC[n, c, lp], zVar], #] & /@ s, dkc = Sum[gt[l, lp, L] daCv[[L + 1]], {L, Range[Abs[l - lp], l + lp]}],
          tc = Sum[gt[l, lp, L] angF[l, lp, L] aCv[[L + 1]], {L, Range[Abs[l - lp], l + lp]}]},
         -2 (Im[dkc] Im[drcv] + Im[tc] Im[rcv]/s^2)]];
      {Total[w g0 fb], Total[w g1 fb]}], {lp, 0, lpMax}]];
homOnGrid[l_, aRv_, {u_, du_}, grid_] := Module[{AB = LinearSolve[{{aRv^l, aRv^(-l - 1)}, {l aRv^(l - 1), -(l + 1) aRv^(-l - 2)}}, {u, du}]},
   (AB[[1]] #^l + AB[[2]] #^(-l - 1)) & /@ grid];

(* ================= main loop ================= *)
rows = {};      (* {case, n, Rscheme, l'max, l, total, b, c} scaled *)
Do[
  SetCSParameters[casePars[Q]];
  Do[Module[{c = r0^-(n + 1), nodes, mv, scale, exactSrc, cCurves = <||>, bFixed = <||>},
     nodes = cNodes[12];
     mv = modeVals[Q, n, c, #] & /@ nodes[[All, 1]];
     Print["    ", Q, " n=", n, ": node values cached [", tick[], " s]"];
     Do[
      scale[l] = Max[Abs[analyticEval[totalSol[Q, l], #] & /@ gCom]];
      exactSrc = Table[Module[{s = nodes[[k, 1]]}, pwVal[Sum[pieceAnalyticSource[Q, p, l], {p, {"SS", "SR", "RR"}}], s] - pwVal[{ssPsrc[Q, n, c, l, "real"], ssPsrc[Q, n, c, l, "real"]}, s]], {k, Length[nodes]}];
      Do[cCurves[{l, M}] = cResponse[l, nodes, Table[truncSrc[Q, l, M, mv[[k]], nodes[[k, 1]]], {k, Length[nodes]}] - exactSrc, gCom], {M, lpList}];
      Do[bFixed[{l, R}] = Accumulate[bPerLp[Q, n, c, R, l, lpTop, 30]], {R, Rfixed}];
      Do[
       Module[{cc = cCurves[{l, M}], schemes},
        schemes = Join[Table[{ToString[R, InputForm], R, bFixed[{l, R}][[M + 1]]}, {R, Rfixed}],
          {{"r0/l'max", r0/M, Total[bPerLp[Q, n, c, r0/M, l, M, 30]]}}];
        Do[Module[{bh = Join[ConstantArray[0, Length[gIn]], homOnGrid[l, r0 + sc[[2]], sc[[3]], gOut]]},
           AppendTo[rows, {Q, n, sc[[1]], M, l, Max[Abs[cc + bh]]/scale[l], Max[Abs[bh]]/scale[l], Max[Abs[cc]]/scale[l]}]], {sc, schemes}]],
       {M, lpList}],
      {l, lOut}];
     Print["    ", Q, " n=", n, " done [", tick[], " s]"]], {n, nS}],
  {Q, {"Product", "GradDot"}}];

(* ================= verification against Stage 6 machinery ================= *)
Print["\n--- verification ---"];
Do[SetCSParameters[casePars[Q]];
  Module[{n = 3, c = r0^-4, M = 8, nodes = cNodes[12], nodes2 = cNodes[20], l = 3},
   (* V0: quadrature convergence of (c) *)
   Module[{d1, d2, ex1, ex2},
    ex1 = Table[truncSrc[Q, l, M, modeVals[Q, n, c, s], s] - pwVal[Sum[pieceAnalyticSource[Q, p, l], {p, {"SS", "SR", "RR"}}], s] + pwVal[{ssPsrc[Q, n, c, l, "real"], ssPsrc[Q, n, c, l, "real"]}, s], {s, nodes[[All, 1]]}];
    ex2 = Table[truncSrc[Q, l, M, modeVals[Q, n, c, s], s] - pwVal[Sum[pieceAnalyticSource[Q, p, l], {p, {"SS", "SR", "RR"}}], s] + pwVal[{ssPsrc[Q, n, c, l, "real"], ssPsrc[Q, n, c, l, "real"]}, s], {s, nodes2[[All, 1]]}];
    d1 = cResponse[l, nodes, ex1, gCom]; d2 = cResponse[l, nodes2, ex2, gCom];
    csCheck[Q <> ": (c) quadrature converged (12 vs 20 nodes per subinterval, relative)", Max[Abs[d1 - d2]]/Max[Abs[d2]], 10^-10]];
   (* V1: fast (b) = Stage 6 decomp (same nodes), R = 1/2 *)
   Block[{aL = r0 - 1/2, aR = r0 + 1/2},
    Module[{fast = bPerLp[Q, n, c, 1/2, l, 12, 30], ref = Table[decomp[Q, n, c, l, {lp}, Join[gl[aL, r0, 30], gl[r0, aR, 30]]]["b"], {lp, {0, 5, 12}}]},
     csCheck[Q <> ": fast (b) per l' equals Stage 6 decomp (l' = 0, 5, 12)", Max[Abs[fast[[{1, 6, 13}]] - ref]]/Max[Abs[ref]], 10^-12]]];
   (* V2 and V3: full Stage 6 runs (three-region contour + real-line solves) vs decomposition, on the Stage 6 grid of each tube *)
   Do[Module[{R = vr[[1]], Mv = vr[[2]], res, mine, sc},
      SetCSParameters[Join[casePars[Q], <|"radius" -> R, "amp" -> 2 R|>]];
      Block[{aL = r0 - R, aR = r0 + R, extPts = Join[Subdivide[rMin, r0 - R, 35], Subdivide[r0 + R, rMax, 45]], extR = Subdivide[r0 + R, rMax, 45]},
       res = runSetting[Q, n, c, Mv, {l}];
       Module[{g = Join[Subdivide[rMin, r0 - R, 35], Subdivide[r0 + R, rMax, 45]], nd, ex, cc, bh},
        nd = cNodes[12, g];
        ex = Table[truncSrc[Q, l, Mv, modeVals[Q, n, c, s], s] - pwVal[Sum[pieceAnalyticSource[Q, p, l], {p, {"SS", "SR", "RR"}}], s] + pwVal[{ssPsrc[Q, n, c, l, "real"], ssPsrc[Q, n, c, l, "real"]}, s], {s, nd[[All, 1]]}];
        cc = cResponse[l, nd, ex, g];
        bh = Join[ConstantArray[0, 36], homOnGrid[l, r0 + R, Total[bPerLp[Q, n, c, R, l, Mv, 30]], Subdivide[r0 + R, rMax, 45]]];
        sc = res[l]["scale"];
        csCheck[Q <> ": R=" <> ToString[R, InputForm] <> ", l'max=" <> ToString[Mv] <> ": (c) by quadrature = Stage 6 real-line (c) (max over grid, relative)", Abs[Max[Abs[cc]]/sc - res[l]["Ec"]]/res[l]["Ec"], 10^-6];
        csCheck[Q <> ": R=" <> ToString[R, InputForm] <> ", l'max=" <> ToString[Mv] <> ": (b)+(c) = total from the contour run (max over grid, relative)", Abs[Max[Abs[cc + bh]]/sc - res[l]["Etot"]]/res[l]["Etot"], 10^-6]]]],
    {vr, If[smoke, {{1/2, 8}}, {{1/2, 8}, {r0/12, 12}}]}];
   SetCSParameters[casePars[Q]]], {Q, {"Product", "GradDot"}}];

(* ================= tables, figure ================= *)
mxr[Q_, n_, sch_, M_, k_] := Max[Cases[rows, {Q, n, sch, M, __}][[All, k]]];
schemes = {"1/2", "1/4", "r0/l'max"};
md = StringRiffle[Join[{"| case | n | R | quantity | " <> StringRiffle[ToString /@ lpList, " | "] <> " |", "|---|---|---|---|" <> StringRepeat["---|", Length[lpList]]},
    Flatten@Table[Table["| " <> Q <> " | " <> ToString[n] <> " | " <> sch <> " | " <> {"total", "(b)", "(c)"}[[k - 5]] <> " | " <> StringRiffle[sci[mxr[Q, n, sch, #, k]] & /@ lpList, " | "] <> " |",
       {k, 6, 8}], {Q, {"Product", "GradDot"}}, {n, nS}, {sch, schemes}]], "\n"];
optRows = Flatten[Table[Module[{tot = mxr[Q, n, sch, #, 6] & /@ lpList, k},
     k = First[Ordering[tot, 1]]; {Q, n, sch, lpList[[k]], tot[[k]], tot[[-1]], If[OrderedQ[Reverse[tot]], "monotone decreasing", "not monotone"]}], {Q, {"Product", "GradDot"}}, {n, nS}, {sch, schemes}], 2];
md2 = StringRiffle[Join[{"| case | n | R | argmin l'max | min total | total at l'max = 40 | behaviour |", "|---|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> #[[3]] <> " | " <> ToString[#[[4]]] <> " | " <> sci[#[[5]]] <> " | " <> sci[#[[6]]] <> " | " <> #[[7]] <> " |") & /@ optRows], "\n"];
Print["\n", md2];
If[smoke, csSummary[]; Exit[]];
Export[FileNameJoin[{resDir, "stage7_summary.md"}],
  "# Stage 7: l'max optimum and shrinking worldtube\n\nScaled exterior error (max over l = 0..6 of max|error|/max|analytic total| over r in [1,3.5] U [6.5,10]); total = (b) + (c) from the validated Lemma-1 decomposition.\n\n" <>
   md2 <> "\n\n## Errors against l'max\n\n" <> md <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage7_errors.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]]} & /@ rows, {"case", "n", "R", "lpmax", "l", "total", "b", "c"}], "CSV"];
$fs = {FontFamily -> "Times", FontSize -> 9};
schStyle = <|"1/2" -> RGBColor[0.1, 0.35, 0.75], "1/4" -> RGBColor[0.85, 0.45, 0.1], "r0/l'max" -> RGBColor[0.15, 0.6, 0.25]|>;
panels = Table[ListLogPlot[Join[Flatten[Table[Table[{M, mxr[Q, n, sch, M, 6]}, {M, lpList}], {n, nS}, {sch, schemes}], 1], Table[Table[{M, mxr[Q, n, "1/2", M, 8]}, {M, lpList}], {n, nS}]],
    Joined -> True, PlotMarkers -> Join[Flatten[Table[{{"\[FilledCircle]", 5}, {"\[FilledSquare]", 5}, {"\[FilledDiamond]", 6}}, {n, nS}], 1], {None, None}],
    PlotStyle -> Join[Flatten[Table[Directive[schStyle[sch], AbsoluteThickness[1.3], If[n == 3, Dashing[{}], Dashed]], {n, nS}, {sch, schemes}], 1],
      Table[Directive[Gray, AbsoluteThickness[0.8], If[n == 3, Dotted, DotDashed]], {n, nS}]],
    Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]'\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 212,
    Epilog -> {Inset[Style[If[Q === "Product", "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(2\)]\)", "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(4\)]\)"], $fs], Scaled[{0.97, 0.97}], {Right, Top}]}], {Q, {"Product", "GradDot"}}];
leg = Column[{LineLegend[Values[schStyle], {"R = 1/2", "R = 1/4", "R = \!\(\*SubscriptBox[\(r\), \(0\)]\)/\[ScriptL]'\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)"}, LegendLayout -> "Row", LabelStyle -> $fs],
    LineLegend[{Directive[Black], Directive[Black, Dashed], Directive[Gray, Dotted], Directive[Gray, DotDashed]}, {"total, n = 3", "total, n = 5", "(c) only, n = 3", "(c) only, n = 5"}, LegendLayout -> "Row", LabelStyle -> $fs]}, Alignment -> Center];
Export[FileNameJoin[{figDir, "fig_lprime_optimum_cauchy.pdf"}], Legended[Grid[{panels}, Spacings -> 1], Placed[leg, Below]]];
Export[FileNameJoin[{figDir, "fig_lprime_optimum_cauchy.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], N[#[[6]], 6], N[#[[7]], 6], N[#[[8]], 6]} & /@ rows, {"case", "n", "R", "lpmax", "l", "total", "b", "c"}], "CSV"];
Print["\nStage 7 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage7_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
