#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Figures for the LaTeX note (written to ../figures/, data as CSV next to each PDF).
   Both cases: Q = "Product" (1/rho^2; PhiR = c0 + c1 r cos, Stage 3 parameters) and Q = "GradDot" (1/rho^4; plus c2 r^2 P2,
   Stage 5 parameters).  K continued by Theta -> H, Cauchy data at r = 1 from the analytic particular solutions.
   Run:  wolframscript -file make_figures.wls [lmaxTot]      (default 20)                                              *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
figDir = FileNameJoin[{ParentDirectory[dir], "figures"}]; Quiet[CreateDirectory[figDir]];
lmaxTot = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], 20];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
cases = <|"Product" -> <|"pars" -> <|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>, "lps" -> {0, 1}|>,
   "GradDot" -> <|"pars" -> <|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10, "c2" -> 1/100|>, "lps" -> {0, 1, 2}|>|>;
cases["Product", "label"] = "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(2\)]\)  (Q[A,B] = AB)";
cases["GradDot", "label"] = "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(4\)]\)  (Q[A,B] = \[Del]A\[CenterDot]\[Del]B)";
pieces = {"SS", "SR", "RR"};

(* ---------------- compute all modes for one case ---------------- *)
computeCase[Q_] := Module[{c = cases[Q], aL, aR, rMin, rMax, res1, ext, til, line, out = <||>},
   SetCSParameters[c["pars"]];
   aL = tubeL[]; aR = tubeR[]; rMin = csPar["rMin"]; rMax = csPar["rMax"];
   res1 = Association@Table[lp -> firstOrderRegular[lp], {lp, c["lps"]}];
   til = Association@Table[lp -> res1[lp]["sol2"], {lp, c["lps"]}];
   If[Q === "Product",
    ext = Association@Table[With[{lp = lp}, lp -> Module[{f}, f[r_?NumericQ] := exteriorValue[res1[lp], r]; f]], {lp, c["lps"]}];
    line = Association@Table[lp -> realLineSolve[lp, Function[r, 0], {phiRAnalytic[lp][rMin], phiRAnalytic[lp]'[rMin]}], {lp, c["lps"]}],
    ext = Association@Table[With[{lp = lp}, lp -> {Module[{f}, f[r_?NumericQ] := exteriorValue[res1[lp], r]; f],
          Module[{g}, g[r_?NumericQ] := If[r <= aL, res1[lp]["sol1"]'[r], res1[lp]["sol3"]'[r]]; g]}], {lp, c["lps"]}];
    line = Association@Table[With[{s = realLineSolve[lp, Function[r, 0], {phiRAnalytic[lp][rMin], phiRAnalytic[lp]'[rMin]}]}, lp -> {s, s'}], {lp, c["lps"]}]];
   Do[
    Module[{sol = Association@Table[p -> pieceAnalyticSolution[Q, p, l], {p, pieces}], d, rS, rR, lc},
     d[p_] := {analyticEval[sol[p], rMin], analyticEvalD[sol[p], rMin]};
     rS = threeRegionSolve[l, ssSourceReal[Q, l], ssSourceContour[Q, l], d["SS"]];
     rR = threeRegionSolve[l, srSourceReal[Q, l, ext], srSourceContour[Q, l, til, "ThetaToH"], d["SR"]];
     lc = realLineSolve[l, rrSourceReal[Q, l, line], d["RR"]];
     out[l] = <|"num" -> Function[r, exteriorValue[rS, r] + exteriorValue[rR, r] + lc[r]],
        "ana" -> Function[r, Total[analyticEval[sol[#], r] & /@ pieces]],
        "numPiece" -> <|"SS" -> Function[r, exteriorValue[rS, r]], "SR" -> Function[r, exteriorValue[rR, r]], "RR" -> Function[r, lc[r]]|>,
        "anaPiece" -> AssociationMap[Function[p, Function[r, analyticEval[sol[p], r]]], pieces]|>;
     If[Mod[l, 4] == 0, Print["    ", Q, " l = ", l, "  [", tick[], " s]"]]],
    {l, 0, lmaxTot}];
   out];

exact3D[Q_, r_, th_] := If[Q === "Product", Total[analytic3D[#, r, th] & /@ pieces], Total[analytic3D["GradDot", #, r, th] & /@ pieces]];
Ylm[l_, th_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, Cos[th]];

(* ---------------- style ---------------- *)
$fs = {FontFamily -> "Times", FontSize -> 9};
cols = Table[ColorData[{"Rainbow", {0, 1}}][x], {x, Subdivide[0, 1, 6]}];
tubeShade[aL_, aR_] := {GrayLevel[0.9], EdgeForm[], Rectangle[{aL, -10^6}, {aR, 10^6}]};
panelW = 212;

data = <||>;
Do[
  Print["computing ", Q, " (l = 0..", lmaxTot, ")"];
  data[Q] = computeCase[Q], {Q, Keys[cases]}];
Print["modes computed [", tick[], " s]"];

(* ================= Figure 1: recovered vs analytic total modes ================= *)
fig1L = {0, 1, 2, 4};
csv1 = {{"case", "l", "r", "numerical_exterior", "analytic"}};
fig1Panels = Table[
   SetCSParameters[cases[Q, "pars"]];
   Module[{aL = tubeL[], aR = tubeR[], r0 = csPar["r0"], gExt, gAll, curvesN, curvesA, yr, legs},
    gExt = {Subdivide[1, aL, 60], Subdivide[aR, 10, 80]};
    gAll = Select[Subdivide[1, 10, 360], Abs[# - r0] > 1/200 &];
    curvesN = Table[Table[{r, data[Q][l]["num"][r]}, {r, seg}], {l, fig1L}, {seg, gExt}];
    curvesA = Table[{r, data[Q][l]["ana"][r]}, {l, fig1L}, {r, gAll}];
    Do[Do[AppendTo[csv1, {Q, fig1L[[k]], N[p[[1]], 8], N[p[[2]], 16], N[data[Q][fig1L[[k]]]["ana"][p[[1]]], 16]}], {seg, curvesN[[k]]}, {p, seg}], {k, Length[fig1L]}];
    yr = MinMax[Flatten[curvesN[[All, All, All, 2]]]]; yr = yr + {-1, 1} 0.08 (yr[[2]] - yr[[1]]);
    Show[
     ListLinePlot[Flatten[curvesN, 1], PlotStyle -> Flatten[Table[{Directive[cols[[2 k - 1]], AbsoluteThickness[1.6]], Directive[cols[[2 k - 1]], AbsoluteThickness[1.6]]}, {k, Length[fig1L]}], 1]],
     ListLinePlot[curvesA, PlotStyle -> Table[Directive[Black, Dashed, AbsoluteThickness[0.8]], {Length[fig1L]}]],
     Prolog -> tubeShade[aL, aR], PlotRange -> {{1, 10}, yr}, Frame -> True, Axes -> False,
     FrameLabel -> {"r", "\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(2,\[ScriptL]\)]\)(r)"}, BaseStyle -> $fs, ImageSize -> panelW,
     Epilog -> {Inset[Style[cases[Q, "label"], $fs], Scaled[{0.03, 0.97}], {Left, Top}]},
     PlotRangePadding -> None]],
   {Q, Keys[cases]}];
leg1 = LineLegend[Join[Table[Directive[cols[[2 k - 1]], AbsoluteThickness[1.6]], {k, Length[fig1L]}], {Directive[Black, Dashed]}],
   Join[("\[ScriptL] = " <> ToString[#]) & /@ fig1L, {"analytic"}], LegendLayout -> "Row", LabelStyle -> $fs];
fig1 = Legended[Grid[{fig1Panels}, Spacings -> 1], Placed[leg1, Below]];
Export[FileNameJoin[{figDir, "fig_modes_recovered.pdf"}], fig1];
Export[FileNameJoin[{figDir, "fig_modes_recovered.csv"}], csv1, "CSV"];
Print["fig 1 written [", tick[], " s]"];

(* ================= Figure 2: pointwise relative and scaled errors ================= *)
fig2L = Range[0, Min[12, lmaxTot]];
csv2 = {{"case", "l", "r", "abs_err", "rel_err", "scaled_err"}};
fig2Panels = Table[
   SetCSParameters[cases[Q, "pars"]];
   Module[{aL = tubeL[], aR = tubeR[], gExt, rel, scl, mk},
    gExt = {Subdivide[1, aL, 70], Subdivide[aR, 10, 90]};
    rel = <||>; scl = <||>;
    Do[Module[{num = data[Q][l]["num"], ana = data[Q][l]["ana"], mx},
       mx = Max[Abs[ana /@ Flatten[gExt]]];
       rel[l] = Table[{r, Abs[num[r] - ana[r]]/Abs[ana[r]] + 10^-30}, {seg, gExt}, {r, seg}];
       scl[l] = Table[{r, Abs[num[r] - ana[r]]/mx + 10^-30}, {seg, gExt}, {r, seg}];
       Do[AppendTo[csv2, {Q, l, N[r, 8], N[Abs[num[r] - ana[r]], 6], N[Abs[num[r] - ana[r]]/Abs[ana[r]], 6], N[Abs[num[r] - ana[r]]/mx, 6]}], {seg, gExt}, {r, seg}]], {l, fig2L}];
    mk[set_, ylab_] := ListLogPlot[Flatten[Table[set[l], {l, fig2L}], 1], Joined -> True,
      PlotStyle -> Flatten[Table[{Directive[ColorData[{"Rainbow", {0, Max[fig2L]}}][l], AbsoluteThickness[1]], Directive[ColorData[{"Rainbow", {0, Max[fig2L]}}][l], AbsoluteThickness[1]]}, {l, fig2L}], 1],
      Prolog -> tubeShade[aL, aR], PlotRange -> {{1, 10}, {10^-22, 10^-8}}, Frame -> True, Axes -> False, FrameLabel -> {"r", ylab},
      BaseStyle -> $fs, ImageSize -> panelW, Epilog -> {Inset[Style[cases[Q, "label"], $fs], Scaled[{0.03, 0.97}], {Left, Top}]}, PlotRangePadding -> None];
    {mk[rel, "|\[CapitalDelta]|/|analytic|"], mk[scl, "|\[CapitalDelta]|/\!\(\*SubscriptBox[\(max\), \(r\)]\)|analytic|"]}],
   {Q, Keys[cases]}];
leg2 = BarLegend[{ColorData[{"Rainbow", {0, Max[fig2L]}}], {0, Max[fig2L]}}, LegendLabel -> "\[ScriptL]", LabelStyle -> $fs, LegendMarkerSize -> {8, 120},
   Ticks -> Range[0, Max[fig2L], 2]];
fig2 = Legended[Grid[Transpose[fig2Panels], Spacings -> {1, 1}], Placed[leg2, Right]];
Export[FileNameJoin[{figDir, "fig_relative_error.pdf"}], fig2];
Export[FileNameJoin[{figDir, "fig_relative_error.csv"}], csv2, "CSV"];
Print["fig 2 written [", tick[], " s]"];

(* ================= Figure 3: convergence with the number of modes ================= *)
pts = {2, 3, 4, 6, 7, 8}; th3 = Pi/3; maxDiff = <||>;
csv3 = {{"case", "r", "theta", "lmax", "err_numerical_modes", "err_analytic_modes", "reference_(r</r>)^lmax"}};
fig3Panels = Table[
   SetCSParameters[cases[Q, "pars"]];
   Module[{r0 = csPar["r0"], curves = {}, styles = {}},
    Do[
     Module[{r = pts[[k]], ex, nm, an, eN, eA, ref, ratio},
      ex = N[exact3D[Q, r, th3], 40];
      nm = Accumulate[Table[data[Q][l]["num"][r] N[Ylm[l, th3], 40], {l, 0, lmaxTot}]];
      an = Accumulate[Table[data[Q][l]["ana"][r] N[Ylm[l, th3], 40], {l, 0, lmaxTot}]];
      eN = Abs[nm - ex] + 10^-30; eA = Abs[an - ex] + 10^-30;
      maxDiff[Q] = Max[Lookup[maxDiff, Q, 0], Max[Abs[nm - an]]];
      ratio = Min[r, r0]/Max[r, r0];
      ref = Table[eA[[5]] ratio^(lm - 4), {lm, 0, lmaxTot}];
      Do[AppendTo[csv3, {Q, r, N[th3, 6], lm, N[eN[[lm + 1]], 6], N[eA[[lm + 1]], 6], N[ref[[lm + 1]], 6]}], {lm, 0, lmaxTot}];
      curves = Join[curves, {Transpose[{Range[0, lmaxTot], eN}], Transpose[{Range[0, lmaxTot], eA}], Transpose[{Range[0, lmaxTot], ref}]}];
      styles = Join[styles, {Directive[cols[[k]], AbsoluteThickness[1.4]], Directive[cols[[k]], AbsoluteThickness[0]], Directive[cols[[k]], Dotted, AbsoluteThickness[0.8]]}]],
     {k, Length[pts]}];
    ListLogPlot[curves, Joined -> Flatten[Table[{True, False, True}, {Length[pts]}]], PlotStyle -> styles,
     PlotMarkers -> Flatten[Table[{None, {"\[EmptyCircle]", 6}, None}, {Length[pts]}], 1],
     PlotRange -> {{0, lmaxTot}, All}, Frame -> True, Axes -> False,
     FrameLabel -> {"\!\(\*SubscriptBox[\(\[ScriptL]\), \(max\)]\)", "|partial sum \[Minus] \!\(\*SubsuperscriptBox[\(\[CapitalPhi]\), \(2\), \(exact\)]\)|"},
     BaseStyle -> $fs, ImageSize -> panelW, Epilog -> {Inset[Style[cases[Q, "label"], $fs], Scaled[{0.97, 0.97}], {Right, Top}]}]],
   {Q, Keys[cases]}];
leg3 = Column[{LineLegend[Table[Directive[cols[[k]], AbsoluteThickness[1.4]], {k, Length[pts]}], ("r = " <> ToString[#]) & /@ pts, LegendLayout -> "Row", LabelStyle -> $fs],
    Row[{LineLegend[{Directive[Black, AbsoluteThickness[1.4]]}, {"numerical modes"}, LabelStyle -> $fs],
      PointLegend[{Black}, {"analytic modes"}, LegendMarkers -> {{"\[EmptyCircle]", 8}}, LabelStyle -> $fs],
      LineLegend[{Directive[Black, Dotted]}, {"\[Proportional] (\!\(\*SubscriptBox[\(r\), \(<\)]\)/\!\(\*SubscriptBox[\(r\), \(>\)]\)\!\(\*SuperscriptBox[\()\), SubscriptBox[\(\[ScriptL]\), \(max\)]]\), \[Theta] = \[Pi]/3"}, LabelStyle -> $fs]}, "  "]}, Alignment -> Center];
fig3 = Legended[Grid[{fig3Panels}, Spacings -> 1], Placed[leg3, Below]];
Export[FileNameJoin[{figDir, "fig_mode_convergence.pdf"}], fig3];
Export[FileNameJoin[{figDir, "fig_mode_convergence.csv"}], csv3, "CSV"];
Print["fig 3 written [", tick[], " s]"];

(* short numerical summary for the captions *)
Do[SetCSParameters[cases[Q, "pars"]];
  Print[Q, ": max pointwise scaled error (l <= ", Max[fig2L], ") = ", sci[Max[ToExpression /@ Cases[csv2, {Q, __}][[All, 6]]]],
   ";  max over points and lmax <= ", lmaxTot, " of |numerical - analytic partial sum|: ",
   sci[maxDiff[Q]]], {Q, Keys[cases]}];
Print["done [", tick[], " s]"];
