#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 5b: the 1/rho^4 case, Q[A,B] = Grad A . Grad B.  SS and SR on the contour inside the tube (real line outside),
   RR on the real line; Cauchy and two-point BCs; exterior errors and shadows; SS regularisation test (point-supported
   terms at x0?) with positive controls; contour-amplitude independence; 3D reconstruction.
   Run:  wolframscript -file stage5b_gradient_solve.wls [lmax]      (default lmax = 6)                               *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
base = <|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10, "c2" -> 1/100|>;
SetCSParameters[base];
Q = "GradDot";
lmax = If[Length[$ScriptCommandLine] >= 2, ToExpression[$ScriptCommandLine[[2]]], 6];
tag = If[lmax == 6, "", "_lmax" <> ToString[lmax]];
r0 = csPar["r0"]; q = csPar["q"]; rMin = csPar["rMin"]; rMax = csPar["rMax"]; aL = tubeL[]; aR = tubeR[];
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
extPts = Join[Subdivide[rMin, aL, 35], Subdivide[aR, rMax, 45]];
tolExt = 10^-10;
lps = {0, 1, 2};

Print["=== Stage 5b: 1/rho^4 case (Q = Grad.Grad), l = 0..", lmax, ", PhiR = c0 + c1 r cos + c2 r^2 P2, c = {1, 1/10, 1/100} ==="];

(* ---------------- first-order regular field (Stage 2 pipeline), now with l' = 0, 1, 2 ---------------- *)
res1 = Association@Table[lp -> firstOrderRegular[lp], {lp, lps}];
Do[csCheck["first-order phiR_" <> ToString[lp] <> ": contour solution vs continuation (scaled)",
   scaledError[res1[lp]["sol2"] /@ Subdivide[aL, aR, 50], N[phiRAnalytic[lp][yC[#]], 40] & /@ Subdivide[aL, aR, 50]], 10^-12];
  csCheck["first-order phiR_" <> ToString[lp] <> ": Im at aR (value + derivative)", Abs[res1[lp]["imAR"][[1]]] + Abs[res1[lp]["imAR"][[2]]], 10^-15], {lp, lps}];
phiRExt = Association@Table[With[{lp = lp}, lp -> {Module[{f}, f[r_?NumericQ] := exteriorValue[res1[lp], r]; f],
      Module[{g}, g[r_?NumericQ] := If[r <= aL, res1[lp]["sol1"]'[r], res1[lp]["sol3"]'[r]]; g]}], {lp, lps}];
phiRTilde = Association@Table[lp -> res1[lp]["sol2"], {lp, lps}];
phiRLine = Association@Table[With[{s = realLineSolve[lp, Function[r, 0], {phiRAnalytic[lp][rMin], phiRAnalytic[lp]'[rMin]}]}, lp -> {s, s'}], {lp, lps}];
Print["    first-order inputs ready [", tick[], " s]"];

(* ---------------- helpers ---------------- *)
threeEval[res_, r_] := exteriorValue[res, r];
extErrors[res_, sol_, evalF_] := Module[{num = evalF[res, #] & /@ extPts, ref = analyticEval[sol, #] & /@ extPts},
   {scaledError[num, ref], Max[Abs[num - ref]]}];
rows = {}; store = <||>;
runThree[name_, l_, srcR_, srcC_, sol_] := Module[{d1, vals, rc, rb, ec, eb},
   d1 = {analyticEval[sol, rMin], analyticEvalD[sol, rMin]}; vals = {analyticEval[sol, rMin], analyticEval[sol, rMax]};
   rc = threeRegionSolve[l, srcR, srcC, d1]; ec = extErrors[rc, sol, threeEval];
   rb = threeRegionBVP[l, srcR, srcC, vals]; eb = extErrors[rb, sol, threeEval];
   AppendTo[rows, {name, l, ec[[1]], eb[[1]], rc["imAR"][[1]], rc["imAR"][[2]]}];
   csCheck[name <> " l=" <> ToString[l] <> " exterior, Cauchy data at r=1 (scaled)", ec[[1]], tolExt];
   csCheck[name <> " l=" <> ToString[l] <> " exterior, two-point BVP (scaled)", eb[[1]], tolExt];
   store[{name, l, "C"}] = rc; store[{name, l, "B"}] = rb;
   rc];
(* point-term basis: [1/rho]_l (out) minus its in-branch formula continued (homogeneous), and d/dr0 of it *)
pointBasis[l_] := pointBasis[l] = Module[{rr, R0, b1},
    b1 = Sqrt[4 Pi/(2 l + 1)] (R0^l/rr^(l + 1) - rr^l/R0^(l + 1));
    {Function @@ {{rr}, b1 /. R0 -> r0}, Function @@ {{rr}, D[b1, R0] /. R0 -> r0}}];
fitPoint[l_, d_] := Module[{b = pointBasis[l], M},
   M = {{b[[1]][aR], b[[2]][aR]}, {b[[1]]'[aR], b[[2]]'[aR]}};
   LinearSolve[N[M, 40], d]];
deltaAR[res_, sol_] := {Re[res["valR"]] - analyticEval[sol, aR], Re[res["derR"]] - analyticEvalD[sol, aR]};

(* ---------------- main loop ---------------- *)
regRows = {};
Do[
  Print["\n--- l = ", l, " ---  [", tick[], " s]"];
  Module[{solSS = pieceAnalyticSolution[Q, "SS", l], solSR = pieceAnalyticSolution[Q, "SR", l], solRR = pieceAnalyticSolution[Q, "RR", l],
    rS, rA, rB, dAB, lc, lb, ab, srcRR},
   (* SS *)
   rS = runThree["SS", l, ssSourceReal[Q, l], ssSourceContour[Q, l], solSS];
   ab = fitPoint[l, deltaAR[rS, solSS]];
   AppendTo[regRows, {l, ab[[1]], ab[[2]]}];
   csCheck["SS l=" <> ToString[l] <> " regularisation: no monopole term at x0 (|alpha|)", Abs[ab[[1]]], 10^-12];
   csCheck["SS l=" <> ToString[l] <> " regularisation: no z-dipole term at x0 (|beta|)", Abs[ab[[2]]], 10^-12];
   (* whole complex contour solution = continued (q^2/2)[1/rho^2]_l + homogeneous; shadow at aR predicted *)
   Module[{psiC, dd, hom, aG = Subdivide[aL, aR, 50], pred},
    psiC = Function @@ {{zVar}, canonicalizeLogs[solSS[[2]], zVar]};
    dd = Simplify[ComplexExpand[{analyticEval[solSS, aL] - psiC[aL], analyticEvalD[solSS, aL] - (D[psiC[zVar], zVar] /. zVar -> aL)}]];
    hom = homogeneousWithData[l, aL, dd];
    csCheck["SS l=" <> ToString[l] <> " continued r>r0 formula reproduces the r<r0 formula (Re of matching homogeneous = 0)", Max[Abs[N[Re[dd], 30]]], 10^-25];
    pred[a_] := N[psiC[yC[a]] + hom[yC[a]], 40];
    csCheck["SS l=" <> ToString[l] <> " complex contour solution = continued analytic mode + i*homogeneous (scaled)", scaledError[rS["sol2"] /@ aG, pred /@ aG], 10^-12];
    csCheck["SS l=" <> ToString[l] <> " shadow Im psi(aR) equals the analytic prediction", Abs[rS["imAR"][[1]] - Im[pred[aR]]], 10^-15]];
   (* SR, K, d_rK, K^theta continued by Theta -> H; phiR~ from the contour solve, d_r phiR~ = (d/da)/y' *)
   rA = runThree["SR", l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, phiRTilde, "ThetaToH"], solSR];
   rB = threeRegionSolve[l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, phiRTilde, "CorrectedSingular"],
     {analyticEval[solSR, rMin], analyticEvalD[solSR, rMin]}];
   dAB = scaledError[threeEval[rB, #] & /@ extPts, threeEval[rA, #] & /@ extPts];
   csCheck["SR l=" <> ToString[l] <> " exterior: Theta->H vs corrected-phiS continuation agree (scaled)", dAB, tolExt];
   AppendTo[rows, {"SR(corrected-K)", l, extErrors[rB, solSR, threeEval][[1]], Missing[], rB["imAR"][[1]], rB["imAR"][[2]]}];
   (* RR on the real line *)
   srcRR = rrSourceReal[Q, l, phiRLine];
   lc = realLineSolve[l, srcRR, {analyticEval[solRR, rMin], analyticEvalD[solRR, rMin]}];
   lb = realLineBVP[l, srcRR, {analyticEval[solRR, rMin], analyticEval[solRR, rMax]}];
   AppendTo[rows, {"RR", l, extErrors[lc, solRR, #1[#2] &][[1]], extErrors[lb, solRR, #1[#2] &][[1]], 0, 0}];
   csCheck["RR l=" <> ToString[l] <> " exterior, Cauchy (scaled)", extErrors[lc, solRR, #1[#2] &][[1]], tolExt];
   csCheck["RR l=" <> ToString[l] <> " exterior, BVP (scaled)", extErrors[lb, solRR, #1[#2] &][[1]], tolExt];
   store[{"RR", l, "C"}] = lc; store[{"RR", l, "B"}] = lb;
   Print["    alpha = ", sci[ab[[1]]], ", beta = ", sci[ab[[2]]], ";  Im SS(aR) = ", sci[rS["imAR"][[1]]], ", Im SR(aR) = ", sci[rA["imAR"][[1]]],
    " (Theta->H), ", sci[rB["imAR"][[1]]], " (corrected)"]],
  {l, 0, lmax}];

(* ---------------- positive controls for the regularisation fit ---------------- *)
Print["\n--- Positive controls: add a monopole A delta^3 or a z-dipole A' d_r0 delta^3 at x0 to the SS contour source ---"];
ctlRows = {};
Do[
  Module[{solSS = pieceAnalyticSolution[Q, "SS", l], Yl0 = Sqrt[(2 l + 1)/(4 Pi)], ssc, mono, dip, rM, rD, abM, abD},
   ssc = ssSourceContour[Q, l];
   mono = Module[{f}, f[a_?NumericQ] := ssc[a] + N[Yl0/r0^2 deltaExt[yC[SetPrecision[a, 60]]], 60]; f];       (* A = 1 *)
   dip = Module[{f}, f[a_?NumericQ] := ssc[a] + N[Yl0 (I/(Pi (yC[SetPrecision[a, 60]] - r0)^2)/r0^2 - 2 deltaExt[yC[SetPrecision[a, 60]]]/r0^3), 60]; f];  (* A' = 1: -delta'_ext = +I/(Pi (z-r0)^2) *)
   rM = threeRegionSolve[l, ssSourceReal[Q, l], mono, {analyticEval[solSS, rMin], analyticEvalD[solSS, rMin]}];
   rD = threeRegionSolve[l, ssSourceReal[Q, l], dip, {analyticEval[solSS, rMin], analyticEvalD[solSS, rMin]}];
   abM = fitPoint[l, deltaAR[rM, solSS]]; abD = fitPoint[l, deltaAR[rD, solSS]];
   AppendTo[ctlRows, {l, abM[[1]], abM[[2]], abD[[1]], abD[[2]]}];
   csCheck["control l=" <> ToString[l] <> ": monopole A = 1 recovered as alpha = -1/(4 Pi), beta = 0", Max[Abs[abM - {-1/(4 Pi), 0}]], 10^-12];
   csCheck["control l=" <> ToString[l] <> ": dipole A' = 1 recovered as alpha = 0, beta = -1/(4 Pi)", Max[Abs[abD - {0, -1/(4 Pi)}]], 10^-12]],
  {l, 0, Min[lmax, 6]}];

(* ---------------- contour amplitude ---------------- *)
Print["\n--- SS and SR independent of the contour amplitude (amp = 1/2, 2) ---"];
ampRows = {};
d1Store = Association@Flatten@Table[{{"SS", l} -> deltaAR[store[{"SS", l, "C"}], pieceAnalyticSolution[Q, "SS", l]],
      {"SR", l} -> deltaAR[store[{"SR", l, "C"}], pieceAnalyticSolution[Q, "SR", l]]}, {l, 0, Min[lmax, 6]}];
Do[
  SetCSParameters[Join[base, <|"amp" -> amp|>]];
  Module[{r1 = Association@Table[lp -> firstOrderRegular[lp], {lp, lps}], ext, til},
   ext = Association@Table[With[{lp = lp}, lp -> {Module[{f}, f[r_?NumericQ] := exteriorValue[r1[lp], r]; f],
        Module[{g}, g[r_?NumericQ] := If[r <= aL, r1[lp]["sol1"]'[r], r1[lp]["sol3"]'[r]]; g]}], {lp, lps}];
   til = Association@Table[lp -> r1[lp]["sol2"], {lp, lps}];
   Do[
    Module[{solSS = pieceAnalyticSolution[Q, "SS", l], solSR = pieceAnalyticSolution[Q, "SR", l], rS, rR, dS, dR},
     rS = threeRegionSolve[l, ssSourceReal[Q, l], ssSourceContour[Q, l], {analyticEval[solSS, rMin], analyticEvalD[solSS, rMin]}];
     rR = threeRegionSolve[l, srSourceReal[Q, l, ext], srSourceContour[Q, l, til, "ThetaToH"], {analyticEval[solSR, rMin], analyticEvalD[solSR, rMin]}];
     dS = Max[Abs[deltaAR[rS, solSS] - d1Store[{"SS", l}]]]; dR = Max[Abs[deltaAR[rR, solSR] - d1Store[{"SR", l}]]];
     AppendTo[ampRows, {amp, l, dS, dR, rS["imAR"][[1]], rR["imAR"][[1]]}];
     csCheck["amp=" <> ToString[amp, InputForm] <> " l=" <> ToString[l] <> ": SS exterior data at aR unchanged", dS, 10^-15];
     csCheck["amp=" <> ToString[amp, InputForm] <> " l=" <> ToString[l] <> ": SR exterior data at aR unchanged", dR, 10^-15]],
    {l, 0, Min[lmax, 6]}]];
  Print["    amp = ", amp, " done [", tick[], " s]"], {amp, {1/2, 2}}];
SetCSParameters[base];

(* ---------------- totals and 3D reconstruction ---------------- *)
Print["\n--- Total Phi2 modes and 3D reconstruction (Phi2 = Phi1^2/2 - q PhiR(x0)/rho) ---"];
totMax = Max[Table[Module[{num, ref},
     num = (threeEval[store[{"SS", l, "C"}], #] + threeEval[store[{"SR", l, "C"}], #] + store[{"RR", l, "C"}][#]) & /@ extPts;
     ref = Table[Sum[analyticEval[pieceAnalyticSolution[Q, pc, l], r], {pc, {"SS", "SR", "RR"}}], {r, extPts}];
     scaledError[num, ref]], {l, 0, lmax}]];
csCheck["total SS+SR+RR, all l, exterior (scaled)", totMax, tolExt];
Ylm[l_, th_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, Cos[th]];
recRows = Table[Module[{r = p[[1]], th = p[[2]], num, ana, exact},
    num = Sum[(threeEval[store[{"SS", l, "C"}], r] + threeEval[store[{"SR", l, "C"}], r] + store[{"RR", l, "C"}][r]) N[Ylm[l, th], 30], {l, 0, lmax}];
    ana = Sum[Sum[analyticEval[pieceAnalyticSolution[Q, pc, l], r], {pc, {"SS", "SR", "RR"}}] N[Ylm[l, th], 30], {l, 0, lmax}];
    exact = N[Module[{rho = Sqrt[r^2 + r0^2 - 2 r r0 Cos[th]], pr = csPar["c0"] + csPar["c1"] r Cos[th] + csPar["c2"] r^2 LegendreP[2, Cos[th]]},
       (q/rho + pr)^2/2 - q phiRAtParticle[]/rho], 30];
    csCheck["3D r=" <> ToString[N[r]] <> ": numerical vs analytic mode sum (same lmax, rel)", Abs[(num - ana)/exact], 10^-10];
    {N[r, 4], N[th, 4], N[exact, 12], Abs[num - ana], Abs[ana - exact]}], {p, {{2, Pi/3}, {3, 3/10}, {8, Pi/3}, {19/2, 14/5}}}];

(* ---------------- outputs ---------------- *)
get[name_, l_, k_] := Module[{c = Cases[rows, {name, l, __}]}, If[c === {}, Missing[], c[[1, k]]]];
md = StringRiffle[Join[{"| l | SS Cauchy | SS BVP | SR Cauchy | SR BVP | SR corrected-K | RR Cauchy | RR BVP | Im SS(aR) | Im SR(aR) Theta->H | Im SR(aR) corrected-K | alpha | beta |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|---|"},
    Table["| " <> ToString[l] <> " | " <> StringRiffle[sci /@ {get["SS", l, 3], get["SS", l, 4], get["SR", l, 3], get["SR", l, 4], get["SR(corrected-K)", l, 3],
          get["RR", l, 3], get["RR", l, 4], get["SS", l, 5], get["SR", l, 5], get["SR(corrected-K)", l, 5], regRows[[l + 1, 2]], regRows[[l + 1, 3]]}, " | "] <> " |", {l, 0, lmax}]], "\n"];
Print["\n", md];
Export[FileNameJoin[{resDir, "stage5_exterior_errors" <> tag <> ".md"}],
  "Stage 5 (Q = Grad.Grad). Scaled exterior errors (max|num - analytic|/max|analytic| over [1,4.5] U [5.5,10]), shadows Im psi(aR) (Cauchy runs), and the fitted point-term coefficients alpha (monopole), beta (z-dipole) at x0 for SS.\n\n" <> md <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage5_exterior_errors" <> tag <> ".csv"}],
  Prepend[{#[[1]], #[[2]], sci[#[[3]]], If[MissingQ[#[[4]]], "-", sci[#[[4]]]], sci[#[[5]]], sci[#[[6]]]} & /@ rows, {"piece", "l", "scaled_err_Cauchy", "scaled_err_BVP", "Im_psi(aR)", "Im_psi'(aR)"}], "CSV"];
Export[FileNameJoin[{resDir, "stage5_regularisation" <> tag <> ".csv"}],
  Join[{{"l", "alpha (SS)", "beta (SS)"}}, {#[[1]], sci[#[[2]]], sci[#[[3]]]} & /@ regRows, {{}},
   {{"control l", "alpha (monopole A=1; expect -1/(4Pi))", "beta (monopole)", "alpha (dipole A'=1)", "beta (dipole; expect -1/(4Pi))"}},
   {#[[1]], ToString[NumberForm[#[[2]], 12]], sci[#[[3]]], sci[#[[4]]], ToString[NumberForm[#[[5]], 12]]} & /@ ctlRows], "CSV"];
Export[FileNameJoin[{resDir, "stage5_amplitude" <> tag <> ".csv"}],
  Prepend[{ToString[#[[1]], InputForm], #[[2]], sci[#[[3]]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]]} & /@ ampRows, {"amp", "l", "|delta SS(aR) - amp1|", "|delta SR(aR) - amp1|", "Im SS(aR)", "Im SR(aR)"}], "CSV"];
Export[FileNameJoin[{resDir, "stage5_reconstruction" <> tag <> ".csv"}],
  Prepend[{ToString[#[[1]]], ToString[#[[2]]], ToString[#[[3]]], sci[#[[4]]], sci[#[[5]]]} & /@ recRows, {"r", "theta", "Phi2_exact", "num_minus_analytic_modes", "truncation"}], "CSV"];
errData[p_, l_] := Module[{res = store[{p, l, "C"}], sol = pieceAnalyticSolution[Q, p, l], f},
   f = If[p === "RR", #1[#2] &, threeEval];
   {{#, Abs[f[res, #] - analyticEval[sol, #]] + 10^-28} & /@ Subdivide[rMin, aL, 50], {#, Abs[f[res, #] - analyticEval[sol, #]] + 10^-28} & /@ Subdivide[aR, rMax, 60]}];
Export[FileNameJoin[{resDir, "stage5_exterior_errors" <> tag <> ".pdf"}],
  Row[Table[ListLogPlot[Flatten[Table[errData[p, l], {l, 0, lmax}], 1], Joined -> True, PlotRange -> {10^-28, 10^-6},
     PlotStyle -> Flatten[Table[{ColorData[97][l + 1], ColorData[97][l + 1]}, {l, 0, lmax}], 1],
     PlotLegends -> LineLegend[Table[ColorData[97][l + 1], {l, 0, lmax}], Table["\[ScriptL]=" <> ToString[l], {l, 0, lmax}]],
     PlotLabel -> p <> " (1/\[Rho]^4 case, Cauchy data)", AxesLabel -> {"r", "|error|"}, ImageSize -> 360], {p, {"SS", "SR", "RR"}}]]];
Print["\n   wrote results/stage5_*" <> tag];
Print["\nStage 5b wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage5b_checks" <> tag <> ".csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
