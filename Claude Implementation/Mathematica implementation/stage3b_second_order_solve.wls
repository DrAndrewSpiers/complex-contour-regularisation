#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 3b: second-order pieces SS, SR (contour inside the tube, real line outside) and RR (real line),
   1/rho^2 case, l = 0..lmax.  Exterior errors, shadows at aR, Cauchy-data and two-point BCs, 3D reconstruction.
   Run:  wolframscript -file stage3b_second_order_solve.wls [lmax]      (default lmax = 6)             *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>];
Q = "Product";
lmax = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], 6];
r0 = csPar["r0"]; rMin = csPar["rMin"]; rMax = csPar["rMax"]; aL = tubeL[]; aR = tubeR[];
csResetChecks[];
t0 = AbsoluteTime[];
extL = Subdivide[rMin, aL, 35]; extR = Subdivide[aR, rMax, 45]; extPts = Join[extL, extR];
tolExt = 10^-10;    (* scaled exterior error tolerance, set before the runs from the Stage-2 accuracy (~1e-13) *)
tag = If[lmax == 6, "", "_lmax" <> ToString[lmax]];    (* output-file suffix for non-default lmax *)

Print["=== Stage 3b: second order, Q = ", Q, ", l = 0..", lmax, ", q = ", csPar["q"], ", r0 = ", r0, ", c0 = ", csPar["c0"], ", c1 = ", csPar["c1"], " ==="];
Print["    tube [", aL, ", ", aR, "], contour amp ", csPar["amp"], ", domain [", rMin, ", ", rMax, "], WP ", csPar["wp"],
  ", PrecisionGoal ", csPar["precisionGoal"]];

(* ---------------- first-order regular field (Stage 2 pipeline) ---------------- *)
res1 = Association@Table[lp -> firstOrderRegular[lp], {lp, 0, 1}];
phiRExt = Association@Table[With[{lp = lp}, lp -> Module[{f}, f[r_?NumericQ] := exteriorValue[res1[lp], r]; f]], {lp, 0, 1}];
phiRTilde = Association@Table[lp -> res1[lp]["sol2"], {lp, 0, 1}];
phiRTildeExact = Association@Table[With[{lp = lp}, lp -> Module[{f}, f[a_?NumericQ] := N[phiRAnalytic[lp][yC[SetPrecision[a, 60]]], 50]; f]], {lp, 0, 1}];
phiRLine = Association@Table[lp -> realLineSolve[lp, Function[r, 0], {phiRAnalytic[lp][rMin], phiRAnalytic[lp]'[rMin]}], {lp, 0, 1}];
Do[csCheck["first-order phiR_" <> ToString[lp] <> ": contour solution vs continuation (input to SR)",
   scaledError[phiRTilde[lp] /@ Subdivide[aL, aR, 50], N[phiRAnalytic[lp][yC[#]], 40] & /@ Subdivide[aL, aR, 50]], 10^-12], {lp, 0, 1}];
Print["    first-order inputs ready (", Round[AbsoluteTime[] - t0, 0.1], " s)"];

(* ---------------- helpers ---------------- *)
extErrors[res_, sol_, evalF_] := Module[{num, ref},
   num = evalF[res, #] & /@ extPts; ref = analyticEval[sol, #] & /@ extPts;
   {scaledError[num, ref], Max[Abs[num - ref]], Max[Abs[ref]]}];
threeEval[res_, r_] := exteriorValue[res, r];
lineEval[sol_, r_] := sol[r];

rows = {};          (* {piece, l, bc, scaled, maxAbs, scale, ImVal, ImDer} *)
store = <||>;       (* solutions for plots and the 3D reconstruction *)

runThree[name_, l_, srcR_, srcC_, sol_] := Module[{d1, vals, rc, rb, ec, eb},
   d1 = {analyticEval[sol, rMin], analyticEvalD[sol, rMin]};
   vals = {analyticEval[sol, rMin], analyticEval[sol, rMax]};
   rc = threeRegionSolve[l, srcR, srcC, d1];
   ec = extErrors[rc, sol, threeEval];
   AppendTo[rows, {name, l, "Cauchy@1", ec[[1]], ec[[2]], ec[[3]], rc["imAR"][[1]], rc["imAR"][[2]]}];
   csCheck[name <> " l=" <> ToString[l] <> " exterior, Cauchy data at r=1 (scaled)", ec[[1]], tolExt];
   rb = threeRegionBVP[l, srcR, srcC, vals];
   eb = extErrors[rb, sol, threeEval];
   AppendTo[rows, {name, l, "BVP 1&10", eb[[1]], eb[[2]], eb[[3]], rb["imAR"][[1]], rb["imAR"][[2]]}];
   csCheck[name <> " l=" <> ToString[l] <> " exterior, two-point BVP (scaled)", eb[[1]], tolExt];
   store[{name, l, "C"}] = rc; store[{name, l, "B"}] = rb;
   Print["    ", name, " l=", l, ":  Im psi(aR) = ", sci[rc["imAR"][[1]]], ",  Im psi'(aR) = ", sci[rc["imAR"][[2]]], "   [", Round[AbsoluteTime[] - t0, 0.1], " s]"];
   rc];

(* ---------------- main loop ---------------- *)
Do[
  Print["\n--- l = ", l, " ---"];
  Module[{solSS, solSR, solRR, rS, rA, rB, rX, dAB, dAX, srcRR, lc, lb, ec, eb},
   solSS = pieceAnalyticSolution[Q, "SS", l];
   solSR = pieceAnalyticSolution[Q, "SR", l];
   solRR = pieceAnalyticSolution[Q, "RR", l];

   (* SS: continued source on the contour, real source outside *)
   rS = runThree["SS", l, ssSourceReal[Q, l], ssSourceContour[Q, l], solSS];
   (* Prediction of the whole complex contour solution: the r > r0 formula of [q^2 Log rho]_l with Log[r - r0] -> Log[z - r0]
      (psiC), plus the homogeneous solution that matches the real data at aL.  Its real part must vanish (no D_l-type
      ambiguity), and Im psiC(aR) = 0, so the shadow Im psi(aR) is predicted exactly. *)
   Module[{psiC, dd, hom, aG = Subdivide[aL, aR, 50], pred},
    psiC = Function @@ {{zVar}, canonicalizeLogs[solSS[[2]], zVar]};
    dd = Simplify[ComplexExpand[{analyticEval[solSS, aL] - psiC[aL], analyticEvalD[solSS, aL] - (D[psiC[zVar], zVar] /. zVar -> aL)}]];
    hom = homogeneousWithData[l, aL, dd];
    csCheck["SS l=" <> ToString[l] <> " continued r>r0 formula reproduces the r<r0 formula (Re of matching homogeneous = 0)",
     Max[Abs[N[Re[dd], 30]]], 10^-25];
    pred[a_] := N[psiC[yC[a]] + hom[yC[a]], 40];
    csCheck["SS l=" <> ToString[l] <> " complex contour solution = continued analytic mode + i*homogeneous (scaled)",
     scaledError[rS["sol2"] /@ aG, pred /@ aG], 10^-12];
    csCheck["SS l=" <> ToString[l] <> " shadow Im psi(aR) equals the analytic prediction", Abs[rS["imAR"][[1]] - Im[pred[aR]]], 10^-15]];

   (* SR, method A: K continued by Theta -> H, phiR~ = Stage-2 contour solution *)
   rA = runThree["SR", l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, phiRTilde, "ThetaToH"], solSR];

   (* SR cross-check B: K continued with the corrected singular modes of Stage 1 *)
   rB = threeRegionSolve[l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, phiRTilde, "CorrectedSingular"],
     {analyticEval[solSR, rMin], analyticEvalD[solSR, rMin]}];
   dAB = scaledError[threeEval[rB, #] & /@ extPts, threeEval[rA, #] & /@ extPts];
   AppendTo[rows, {"SR(corrected-K)", l, "Cauchy@1", Sequence @@ extErrors[rB, solSR, threeEval], rB["imAR"][[1]], rB["imAR"][[2]]}];
   csCheck["SR l=" <> ToString[l] <> " exterior: Theta->H vs corrected-singular continuation of K agree (scaled)", dAB, tolExt];
   Print["    SR l=", l, ": shadow at aR, Theta->H: ", sci[rA["imAR"][[1]]], ";  corrected-K: ", sci[rB["imAR"][[1]]], "  (differ, as expected)"];

   (* SR cross-check X: exact analytic continuation of phiR instead of the contour solution *)
   rX = threeRegionSolve[l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, phiRTildeExact, "ThetaToH"],
     {analyticEval[solSR, rMin], analyticEvalD[solSR, rMin]}];
   dAX = scaledError[threeEval[rX, #] & /@ extPts, threeEval[rA, #] & /@ extPts];
   csCheck["SR l=" <> ToString[l] <> " exterior: contour phiR~ vs exact continuation of phiR (scaled)", dAX, tolExt];
   csCheck["SR l=" <> ToString[l] <> " shadow at aR: contour phiR~ vs exact continuation (abs)", Abs[rA["imAR"][[1]] - rX["imAR"][[1]]], 10^-12];

   (* RR: real line throughout, Gaunt coupling of the real-line phiR modes *)
   srcRR = rrSourceReal[Q, l, phiRLine];
   lc = realLineSolve[l, srcRR, {analyticEval[solRR, rMin], analyticEvalD[solRR, rMin]}];
   lb = realLineBVP[l, srcRR, {analyticEval[solRR, rMin], analyticEval[solRR, rMax]}];
   ec = extErrors[lc, solRR, lineEval]; eb = extErrors[lb, solRR, lineEval];
   AppendTo[rows, {"RR", l, "Cauchy@1", ec[[1]], ec[[2]], ec[[3]], 0, 0}];
   AppendTo[rows, {"RR", l, "BVP 1&10", eb[[1]], eb[[2]], eb[[3]], 0, 0}];
   csCheck["RR l=" <> ToString[l] <> " exterior, Cauchy data at r=1 (scaled)", ec[[1]], tolExt];
   csCheck["RR l=" <> ToString[l] <> " exterior, two-point BVP (scaled)", eb[[1]], tolExt];
   csCheck["RR l=" <> ToString[l] <> " inside the tube on the real line (scaled)",
    scaledError[lc /@ Subdivide[aL, aR, 20], analyticEval[solRR, #] & /@ Subdivide[aL, aR, 20]], tolExt];
   store[{"RR", l, "C"}] = lc; store[{"RR", l, "B"}] = lb;
   ], {l, 0, lmax}];

(* ---------------- totals per l ---------------- *)
Print["\n--- Total Phi2 modes (SS + SR + RR), exterior ---"];
totRows = {};
Do[
  Do[
   Module[{num, ref},
    num = (threeEval[store[{"SS", l, bc}], #] + threeEval[store[{"SR", l, bc}], #] + store[{"RR", l, bc}][#]) & /@ extPts;
    ref = (analyticEval[pieceAnalyticSolution[Q, "SS", l], #] + analyticEval[pieceAnalyticSolution[Q, "SR", l], #] + analyticEval[pieceAnalyticSolution[Q, "RR", l], #]) & /@ extPts;
    AppendTo[totRows, {l, bc, scaledError[num, ref], Max[Abs[num - ref]]}];
    csCheck["total l=" <> ToString[l] <> If[bc == "C", " Cauchy", " BVP"] <> " exterior (scaled)", scaledError[num, ref], tolExt]], {bc, {"C", "B"}}], {l, 0, lmax}];

(* ---------------- 3D reconstruction ---------------- *)
Print["\n--- 3D reconstruction at exterior points, Sum_{l<=", lmax, "} phi_l(r) Y_l(theta) ---"];
Ylm[l_, th_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, Cos[th]];
pts3D = {{2, Pi/3}, {3, 3/10}, {4, Pi/2}, {6, 2}, {8, Pi/3}, {19/2, 14/5}};
recRows = Table[
   Module[{r = p[[1]], th = p[[2]], num, numB, ana, exact},
    num = Sum[(threeEval[store[{"SS", l, "C"}], r] + threeEval[store[{"SR", l, "C"}], r] + store[{"RR", l, "C"}][r]) N[Ylm[l, th], 30], {l, 0, lmax}];
    numB = Sum[(threeEval[store[{"SS", l, "B"}], r] + threeEval[store[{"SR", l, "B"}], r] + store[{"RR", l, "B"}][r]) N[Ylm[l, th], 30], {l, 0, lmax}];
    ana = Sum[(analyticEval[pieceAnalyticSolution[Q, "SS", l], r] + analyticEval[pieceAnalyticSolution[Q, "SR", l], r] + analyticEval[pieceAnalyticSolution[Q, "RR", l], r]) N[Ylm[l, th], 30], {l, 0, lmax}];
    exact = N[analytic3D["SS", r, th] + analytic3D["SR", r, th] + analytic3D["RR", r, th], 30];
    csCheck["3D r=" <> ToString[N[r]] <> ", theta=" <> ToString[N[th, 3]] <> ": numerical vs analytic mode sum (same lmax, rel)", Abs[(num - ana)/exact], 10^-10];
    {N[r, 4], N[th, 4], N[exact, 12], N[num, 12], Abs[num - ana], Abs[numB - ana], Abs[ana - exact]}], {p, pts3D}];
Print["    r | theta | Phi2 exact (3D) | mode sum (numerical) | |num - analytic modes| Cauchy | same, BVP | truncation |analytic modes - exact|"];
Scan[Print["    ", StringRiffle[{ToString[#[[1]]], ToString[#[[2]]], ToString[#[[3]]], ToString[#[[4]]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]]}, " | "]] &, recRows];

(* ---------------- negative controls ---------------- *)
Print["\n--- Negative controls: deliberate mistakes that the exterior test must detect ---"];
ctlRows = {};
(* A: SS source continued with Log[(z - r0)^2] instead of 2 Log[z - r0]: 2 Pi i jump on the contour at a = r0 *)
ssWrong[l_] := Module[{zz}, contourFunction[csPar["q"]^2 Sqrt[Pi (2 l + 1)]/(zz r0) legendreQNeumann[l, xOfZ[zz], 2 Log[zz + r0] - Log[(zz - r0)^2]], zz]];
Do[Module[{sol = pieceAnalyticSolution[Q, "SS", l], res, e},
   res = threeRegionSolve[l, ssSourceReal[Q, l], ssWrong[l], {analyticEval[sol, rMin], analyticEvalD[sol, rMin]}];
   e = extErrors[res, sol, threeEval][[1]];
   Print["    SS l=", l, " with Log[(z-r0)^2] on the contour: scaled exterior error = ", sci[e]];
   AppendTo[ctlRows, {"SS, Log[(z-r0)^2] continuation", l, e}];
   csCheck["control A l=" <> ToString[l] <> ": wrong-branch SS source is detected (error > 1e-8)", If[e > 10^-8, 0, 1], 0]], {l, 0, Min[2, lmax]}];
(* B: phiR~_0 from the Theta->H puncture (violates C2: O(1) shadow) used in SR *)
res0bad = firstOrderRegular[0, Automatic, "ThetaToH"];
csCheck["control B: first-order phiR_0 exterior still exact with the Theta->H puncture (Lemma 1)",
  scaledError[exteriorValue[res0bad, #] & /@ extPts, N[phiRAnalytic[0][#], 40] & /@ extPts], 10^-12];
Print["    Theta->H puncture: max |Im phiR~_0| on the contour = ", sci[Max[Abs[Im[res0bad["sol2"] /@ Subdivide[aL, aR, 50]]]]],
  ",  Im phiR~_0(aR) = ", sci[Im[res0bad["sol2"][aR]]]];
(* Detection threshold 1e-8 = 100 x tolExt, i.e. far above the clean-run errors (<= 1e-12).
   (A first version used 1e-3; the l = 0, 1 errors, 2e-5 and 3e-4, fell below it.  See README.) *)
Do[Do[Module[{sol = pieceAnalyticSolution[Q, "SR", l], res, e},
   res = threeRegionSolve[l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, <|0 -> res0bad["sol2"], 1 -> phiRTilde[1]|>, method],
     {analyticEval[sol, rMin], analyticEvalD[sol, rMin]}];
   e = extErrors[res, sol, threeEval][[1]];
   Print["    SR l=", l, " (K continued by ", method, ") with the shadowed phiR~_0: scaled exterior error = ", sci[e]];
   AppendTo[ctlRows, {"SR, phiR~_0 with O(1) shadow (Theta->H puncture), K: " <> method, l, e}];
   csCheck["control B l=" <> ToString[l] <> ", K " <> method <> ": shadowed phiR~ in SR is detected (error > 1e-8)", If[e > 10^-8, 0, 1], 0]],
  {method, {"ThetaToH", "CorrectedSingular"}}], {l, 0, Min[2, lmax]}];

(* ---------------- tables ---------------- *)
Export[FileNameJoin[{resDir, "stage3_exterior_errors" <> tag <> ".csv"}],
  Prepend[{#[[1]], #[[2]], #[[3]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]]} & /@ rows,
   {"piece", "l", "BC", "scaled_ext_err", "max_abs_ext_err", "max|analytic|", "Im_psi(aR)", "Im_psi'(aR)"}], "CSV"];
Export[FileNameJoin[{resDir, "stage3_negative_controls" <> tag <> ".csv"}], Prepend[{#[[1]], #[[2]], sci[#[[3]]]} & /@ ctlRows, {"control", "l", "scaled_ext_err"}], "CSV"];
Export[FileNameJoin[{resDir, "stage3_totals" <> tag <> ".csv"}], Prepend[{#[[1]], #[[2]], sci[#[[3]]], sci[#[[4]]]} & /@ totRows, {"l", "BC", "scaled_err", "max_abs_err"}], "CSV"];
Export[FileNameJoin[{resDir, "stage3_reconstruction" <> tag <> ".csv"}],
  Prepend[{ToString[#[[1]]], ToString[#[[2]]], ToString[#[[3]]], ToString[#[[4]]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]]} & /@ recRows,
   {"r", "theta", "Phi2_exact", "Phi2_modesum_numerical", "num_minus_analytic_modes_Cauchy", "num_minus_analytic_modes_BVP", "truncation_error"}], "CSV"];
(* markdown summary table *)
md = Module[{get, lines},
   get[p_, l_, bc_, k_] := sci[First[Cases[rows, {p, l, bc, s_, a_, sc_, iv_, id_} :> {s, a, sc, iv, id}[[k]]], "-"]];
   lines = Join[{"| l | SS Cauchy | SS BVP | SR Cauchy | SR BVP | SR corrected-K | RR Cauchy | RR BVP | Im SS(aR) | Im SR(aR) Theta->H | Im SR(aR) corrected-K |",
      "|---|---|---|---|---|---|---|---|---|---|---|"},
     Table["| " <> ToString[l] <> " | " <> StringRiffle[{get["SS", l, "Cauchy@1", 1], get["SS", l, "BVP 1&10", 1], get["SR", l, "Cauchy@1", 1], get["SR", l, "BVP 1&10", 1],
          get["SR(corrected-K)", l, "Cauchy@1", 1], get["RR", l, "Cauchy@1", 1], get["RR", l, "BVP 1&10", 1],
          get["SS", l, "Cauchy@1", 4], get["SR", l, "Cauchy@1", 4], get["SR(corrected-K)", l, "Cauchy@1", 4]}, " | "] <> " |", {l, 0, lmax}]];
   StringRiffle[lines, "\n"]];
Export[FileNameJoin[{resDir, "stage3_exterior_errors" <> tag <> ".md"}], "Scaled exterior errors max|num - analytic| / max|analytic| over r in [1, aL] U [aR, 10], and shadows Im psi(aR) (Cauchy-data runs).\n\n" <> md <> "\n", "Text"];
Print["\n", md];

(* ---------------- plots ---------------- *)
errData[p_, l_, bc_] := Module[{res = store[{p, l, bc}], sol = pieceAnalyticSolution[Q, p, l], f},
   f = If[p === "RR", lineEval, threeEval];
   {Select[{#, Abs[f[res, #] - analyticEval[sol, #]] + 10^-28} & /@ Subdivide[rMin, aL, 60], True &],
    {#, Abs[f[res, #] - analyticEval[sol, #]] + 10^-28} & /@ Subdivide[aR, rMax, 80]}];
errPlots = Table[
   ListLogPlot[Flatten[Table[errData[p, l, bc], {l, 0, lmax}], 1], Joined -> True, PlotRange -> {10^-28, 10^-6},
    PlotStyle -> Flatten[Table[{ColorData[97][l + 1], ColorData[97][l + 1]}, {l, 0, lmax}], 1],
    PlotLabel -> p <> If[bc == "C", " (Cauchy data at r = 1)", " (two-point BVP)"], AxesLabel -> {"r", "|error|"}, ImageSize -> 400,
    PlotLegends -> LineLegend[Table[ColorData[97][l + 1], {l, 0, lmax}], Table["\[ScriptL]=" <> ToString[l], {l, 0, lmax}]],
    Epilog -> {Gray, Dashed, Line[{{aL, -64}, {aL, 0}}], Line[{{aR, -64}, {aR, 0}}]}],
   {p, {"SS", "SR", "RR"}}, {bc, {"C", "B"}}];
Export[FileNameJoin[{resDir, "stage3_exterior_errors" <> tag <> ".pdf"}], Grid[errPlots]];
solPlots = Table[
   Module[{sol = pieceAnalyticSolution[Q, p, l], res = store[{p, l, "C"}]},
    Show[Plot[analyticEval[sol, r], {r, rMin, rMax}, PlotStyle -> {Black, Dashed}, PlotRange -> All],
     ListPlot[{{#, threeEval[res, #]} & /@ Subdivide[rMin, aL, 40], {#, threeEval[res, #]} & /@ Subdivide[aR, rMax, 50]}, PlotStyle -> {Red, Red}],
     PlotLabel -> p <> ", \[ScriptL] = " <> ToString[l] <> " (dashed: analytic; dots: contour method)", AxesLabel -> {"r", ""}, ImageSize -> 360]],
   {p, {"SS", "SR"}}, {l, 0, 2}];
Export[FileNameJoin[{resDir, "stage3_exterior_solutions" <> tag <> ".pdf"}], Grid[solPlots]];
srcPlots = Table[
   Module[{f = ssSourceContour[Q, l]},
    Plot[{Re[f[a]], Im[f[a]]}, {a, aL, aR}, PlotStyle -> {Blue, Red}, PlotLegends -> {"Re", "Im"},
     PlotLabel -> "continued SS source on the contour, \[ScriptL] = " <> ToString[l], AxesLabel -> {"a", ""}, ImageSize -> 360]], {l, 0, 2}];
Export[FileNameJoin[{resDir, "stage3_SS_source_on_contour" <> tag <> ".pdf"}], Column[srcPlots]];
shadowPlots = Table[
   Module[{res = store[{p, l, "C"}]},
    Plot[{Re[res["sol2"][a]], Im[res["sol2"][a]]}, {a, aL, aR}, PlotStyle -> {Blue, Red}, PlotLegends -> {"Re", "Im (shadow)"},
     PlotLabel -> p <> " contour solution, \[ScriptL] = " <> ToString[l], AxesLabel -> {"a", ""}, ImageSize -> 360]], {p, {"SS", "SR"}}, {l, 0, 1}];
Export[FileNameJoin[{resDir, "stage3_contour_solutions" <> tag <> ".pdf"}], Grid[shadowPlots]];
Print["\n   wrote results/stage3_*.csv, stage3_exterior_errors.md and stage3_*.pdf"];

Print["\nStage 3b wall time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage3b_checks" <> tag <> ".csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
