#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 6: shadow test with a consistent 3D puncture.
   Puncture PhiP = q/rho + c q rho^n (n odd: non-smooth at x0, mimics a truncation error), c = r0^-(n+1) (and 10x at n = 3).
   Split used consistently (C3): SS = Q[PhiP,PhiP] (closed form), SR = 2Q[PhiP,Phi_res], RR = Q[Phi_res,Phi_res],
   Phi_res = Phi1 - PhiP = PhiR - c q rho^n.  On the contour the residual comes from the first-order contour solve with
   S_eff = delta_ext - L[phiP] = -L[c q cont[rho^n]_l] (Theta -> H continuation), so it carries the truncation-order shadow.
   Both cases: Q = "Product" (1/rho^2) and Q = "GradDot" (1/rho^4).  Error decomposition (Lemma 1):
     (a) real-part mismatch (must be ~0), (b) shadow contamination -2 Sum Im K~P Im phi~res (+ gradient analogue),
     (c) l' truncation of the coupling sums.
   Run:  wolframscript -file stage6_3d_puncture_test.wls                                                              *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
Needs["NumericalDifferentialEquationAnalysis`"];
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
casePars = <|"Product" -> <|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>,
   "GradDot" -> <|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10, "c2" -> 1/100|>|>;
lmax = 6; lpMain = 12; nList = {1, 3, 5, 7};
r0 = 5; q = 1; rMin = 1; rMax = 10;
settings = Join[Table[{n, r0^-(n + 1)}, {n, nList}], {{3, 10 r0^-4}}];
lpScan = {2, 4, 8}; amps = {1/2, 2};
quick = MemberQ[$ScriptCommandLine, "quick"];          (* smoke test: wolframscript -file stage6_3d_puncture_test.wls quick *)
If[quick, lmax = 1; lpMain = 4; nList = {1, 3}; settings = {{1, r0^-2}, {3, r0^-4}, {3, 10 r0^-4}}; lpScan = {2}; amps = {2}];
tag = If[quick, "_quick", ""];
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
SetCSParameters[casePars["Product"]];
aL = tubeL[]; aR = tubeR[];
extPts = Join[Subdivide[rMin, aL, 35], Subdivide[aR, rMax, 45]]; extR = Subdivide[aR, rMax, 45];
Print["=== Stage 6: consistent 3D puncture PhiP = q/rho + c q rho^n, n = ", nList, ", c = r0^-(n+1) (and 10x at n = 3), l = 0..", lmax, ", l'max = ", lpMain, " ==="];

(* ---------------- symbolic building blocks ---------------- *)
rhoN[k_, L_] := rhoN[k, L] = Table[uModeIntegral[L, 0, {"pow", k/2}, zVar, br], {br, {"in", "out"}}];      (* [rho^k]_L, in/out *)
cont[n_, L_] := cont[n, L] = With[{io = rhoN[n, L]}, io[[1]] + (io[[2]] - io[[1]]) HExt[zVar]];           (* Theta -> H *)
dcontPw[n_, L_] := dcontPw[n, L] = With[{io = rhoN[n, L]}, D[io[[1]], zVar] + (D[io[[2]], zVar] - D[io[[1]], zVar]) HExt[zVar]];
gt[l1_, l2_, l3_] := gt[l1, l2, l3] = gaunt[l1, l2, l3];
coulPair[L_] := q Sqrt[4 Pi/(2 L + 1)] {zVar^L/r0^(L + 1), r0^L/zVar^(L + 1)};                                   (* [q/rho]_L in/out *)
(* puncture modes a_L and couplings, real (in/out pairs) and continued *)
aR2[n_, c_, L_] := coulPair[L] + c q rhoN[n, L];
aC[n_, c_, L_] := q Sqrt[4 Pi/(2 L + 1)] coulombRatioContinued[L, zVar, "ThetaToH"] + c q cont[n, L];
daC[n_, c_, L_] := q Sqrt[4 Pi/(2 L + 1)] dCoulombRatioContinued[L, zVar, "ThetaToH"] + c q dcontPw[n, L];
Lrange[l_, lp_] := Range[Abs[l - lp], l + lp];
KPr[n_, c_, l_, lp_] := KPr[n, c, l, lp] = Sum[gt[l, lp, L] aR2[n, c, L], {L, Lrange[l, lp]}];
dKPr[n_, c_, l_, lp_] := dKPr[n, c, l, lp] = D[KPr[n, c, l, lp], zVar];
KthPr[n_, c_, l_, lp_] := KthPr[n, c, l, lp] = Sum[gt[l, lp, L] (L (L + 1) + lp (lp + 1) - l (l + 1))/2 aR2[n, c, L], {L, Lrange[l, lp]}];
KPc[n_, c_, l_, lp_] := KPc[n, c, l, lp] = Sum[gt[l, lp, L] aC[n, c, L], {L, Lrange[l, lp]}];
dKPc[n_, c_, l_, lp_] := dKPc[n, c, l, lp] = Sum[gt[l, lp, L] daC[n, c, L], {L, Lrange[l, lp]}];
KthPc[n_, c_, l_, lp_] := KthPc[n, c, l, lp] = Sum[gt[l, lp, L] (L (L + 1) + lp (lp + 1) - l (l + 1))/2 aC[n, c, L], {L, Lrange[l, lp]}];
(* residual modes: real (in/out) and continued; phiR from the current case parameters *)
resR[n_, c_, lp_] := ConstantArray[phiRAnalytic[lp][zVar], 2] - c q rhoN[n, lp];
resC[n_, c_, lp_] := phiRAnalytic[lp][zVar] - c q cont[n, lp];
pick[{in_, out_}, br_] := If[br === "in", in, out];

(* numeric evaluators *)
evR[e_, r_] := N[e /. zVar -> SetPrecision[r, 70], 50];
pwVal[pair_, r_, br_: Automatic] := evR[pick[pair, If[br === Automatic, If[r < r0, "in", "out"], br]], r];

(* ---------------- SS_P: closed form ---------------- *)
poly[k_, l_] := rhoN[k, l][[1]];                                  (* even k >= 0: polynomial, in = out *)
ssPsrc[Q_, n_, c_, l_, br_] := Which[
   Q === "Product", ssSourceModeExpr[Q, l, zVar, br] + 2 c q^2 poly[n - 1, l] + c^2 q^2 poly[2 n, l],
   True, ssSourceModeExpr[Q, l, zVar, br] - 2 n c q^2 If[n == 1, ssSourceModeExpr["Product", l, zVar, br]/q^2, poly[n - 3, l]] + n^2 c^2 q^2 poly[2 n - 2, l]];
ssPsol[Q_, n_, c_, l_] := Which[
   Q === "Product", pieceAnalyticSolution["Product", "SS", l] + 2 c q^2 poly[n + 1, l]/((n + 1) (n + 2)) + c^2 q^2 poly[2 n + 2, l]/((2 n + 2) (2 n + 3)),
   True, pieceAnalyticSolution["GradDot", "SS", l] - 2 n c q^2 If[n == 1, pieceAnalyticSolution["Product", "SS", l]/q^2, poly[n - 1, l]/((n - 1) n)] +
    n^2 c^2 q^2 poly[2 n, l]/((2 n) (2 n + 1))];
totalSol[Q_, l_] := Sum[pieceAnalyticSolution[Q, p, l], {p, {"SS", "SR", "RR"}}];

(* ---------------- solvers ---------------- *)
lineSolve[l_, srcIn_, srcOut_, {v_, d_}] := Module[{s1, s2},
   s1 = solveRealSegment[l, srcIn, {rMin, r0}, {v, d}];
   s2 = solveRealSegment[l, srcOut, {r0, rMax}, {s1[r0], s1'[r0]}];
   <|"f" -> Function[r, If[r <= r0, s1[r], s2[r]]], "d" -> Function[r, If[r <= r0, s1'[r], s2'[r]]]|>];
resSolve[n_, c_, lp_] := Module[{seff, e = resC[n, c, lp]},
   seff = contourFunction[-c q radialL[lp, cont[n, lp], zVar], zVar];
   solveContourSegment[lp, seff, {e /. zVar -> aL, D[e, zVar] /. zVar -> aL}]];

(* real SR_P and RR_P sources: built once as symbolic expressions per branch (rational functions), then evaluated numerically *)
(* memo key includes c2: phiR (hence the residual) depends on the case parameters *)
dresR[n_, c_, lp_] := dresRm[n, c, lp, csPar["c2"]];
dresRm[n_, c_, lp_, c2_] := dresRm[n, c, lp, c2] = D[#, zVar] & /@ resR[n, c, lp];
srRealExpr[Q_, n_, c_, l_, lps_, br_] := srRealExpr[Q, n, c, l, lps, br] = If[Q === "Product",
    2 Sum[pick[KPr[n, c, l, lp], br] pick[resR[n, c, lp], br], {lp, lps}],
    2 Sum[pick[dKPr[n, c, l, lp], br] pick[dresR[n, c, lp], br] + pick[KthPr[n, c, l, lp], br] pick[resR[n, c, lp], br]/zVar^2, {lp, lps}]];
rrRealExpr[Q_, n_, c_, l_, lps_, br_] := rrRealExpr[Q, n, c, l, lps, br] = Sum[If[gt[l, l1, l2] == 0, 0,
     gt[l, l1, l2] If[Q === "Product", pick[resR[n, c, l1], br] pick[resR[n, c, l2], br],
       pick[dresR[n, c, l1], br] pick[dresR[n, c, l2], br] + (l1 (l1 + 1) + l2 (l2 + 1) - l (l + 1))/2 pick[resR[n, c, l1], br] pick[resR[n, c, l2], br]/zVar^2]],
    {l1, lps}, {l2, lps}];
branchFun[exprF_] := Module[{fi = numericFunction[exprF["in"], zVar], fo = numericFunction[exprF["out"], zVar], f},
   <|"in" -> fi, "out" -> fo, Automatic -> (f[r_?NumericQ] := If[r < r0, fi[r], fo[r]]; f)|>];
srReal[Q_, n_, c_, l_, lps_] := branchFun[srRealExpr[Q, n, c, l, lps, #] &];
rrReal[Q_, n_, c_, l_, lps_] := branchFun[rrRealExpr[Q, n, c, l, lps, #] &];
srContour[Q_, n_, c_, l_, resSols_Association] := Module[{f, kf, kd, kt, yp},
   If[Q === "Product",
    kf = Association@KeyValueMap[#1 -> contourFunction[KPc[n, c, l, #1], zVar] &, resSols];
    f[a_?NumericQ] := 2 Sum[kf[lp][a] resSols[lp][a], {lp, Keys[resSols]}],
    kd = Association@KeyValueMap[#1 -> contourFunction[dKPc[n, c, l, #1], zVar] &, resSols];
    kt = Association@KeyValueMap[#1 -> contourFunction[KthPc[n, c, l, #1]/zVar^2, zVar] &, resSols];
    yp = Module[{a}, numericFunction[yCp[a], a]];
    f[a_?NumericQ] := 2 Sum[kd[lp][a] resSols[lp]'[a]/yp[a] + kt[lp][a] resSols[lp][a], {lp, Keys[resSols]}]];
   f];

(* Lemma-1 decomposition by Gauss-Legendre quadrature on [aL, r0] and [r0, aR]: responses at aR (zero data at aL) *)
nodesFor[np_] := Join[GaussianQuadratureWeights[np, aL, r0, 60], GaussianQuadratureWeights[np, r0, aR, 60]];
decomp[Q_, n_, c_, l_, lps_, nodes_] := Module[{s = nodes[[All, 1]], w = nodes[[All, 2]], g0, g1, fa = 0, fb = 0},
   g0 = (s^2/(2 l + 1)) (aR^l s^(-l - 1) - s^l aR^(-l - 1));
   g1 = (s^2/(2 l + 1)) (l aR^(l - 1) s^(-l - 1) + (l + 1) s^l aR^(-l - 2));
   Do[Module[{rc = evR[resC[n, c, lp], #] & /@ s, rr = pwVal[resR[n, c, lp], #] & /@ s},
      If[Q === "Product",
       Module[{kc = evR[KPc[n, c, l, lp], #] & /@ s, kr = pwVal[KPr[n, c, l, lp], #] & /@ s},
        fa += 2 (Re[kc] Re[rc] - kr rr); fb += -2 Im[kc] Im[rc]],
       Module[{drc = evR[D[resC[n, c, lp], zVar], #] & /@ s, drr = pwVal[dresR[n, c, lp], #] & /@ s,
         dkc = evR[dKPc[n, c, l, lp], #] & /@ s, dkr = pwVal[dKPr[n, c, l, lp], #] & /@ s,
         tc = evR[KthPc[n, c, l, lp], #] & /@ s, tr = pwVal[KthPr[n, c, l, lp], #] & /@ s},
        fa += 2 (Re[dkc] Re[drc] - dkr drr + (Re[tc] Re[rc] - tr rr)/s^2);
        fb += -2 (Im[dkc] Im[drc] + Im[tc] Im[rc]/s^2)]]], {lp, lps}];
   <|"a" -> {Total[w g0 fa], Total[w g1 fa]}, "b" -> {Total[w g0 fb], Total[w g1 fb]}|>];
homErr[l_, {u_, du_}, grid_, scale_] := Module[{AB = LinearSolve[{{aR^l, aR^(-l - 1)}, {l aR^(l - 1), -(l + 1) aR^(-l - 2)}}, {u, du}]},
   Max[Abs[AB[[1]] #^l + AB[[2]] #^(-l - 1)] & /@ grid]/scale];

(* ---------------- one full run: given case, (n, c), l'max ---------------- *)
runSetting[Q_, n_, c_, lpMax_, lList_, opts___] := Module[{lps = Range[0, lpMax], resSols, out = <||>, nodes = nodesFor[60], o = Association[{opts}]},
   resSols = Association@Table[lp -> resSolve[n, c, lp], {lp, lps}];
   If[TrueQ[o["checkResidual"]],
    Do[csCheck[Q <> " n=" <> ToString[n] <> " c=" <> ToString[c, InputForm] <> " l'=" <> ToString[lp] <> ": contour residual = phiR~ - c q cont[rho^n] (scaled)",
       scaledError[resSols[lp] /@ Subdivide[aL, aR, 30], (N[resC[n, c, lp] /. zVar -> #, 40] &) /@ (N[yC[SetPrecision[#, 60]], 60] & /@ Subdivide[aL, aR, 30])], 10^-12],
     {lp, Select[Union[{0, 1, 2, 5, lpMax}], # <= lpMax &]}]];
   Do[
    Module[{tot = totalSol[Q, l], ssE = ssPsol[Q, n, c, l], dSS, dSR, rSS, rSR, lSR, lRR, scale, Etot, Ec, dB, dec, srF, rrF},
     scale = Max[Abs[analyticEval[tot, #] & /@ extPts]];
     dSS = {analyticEval[ssE, rMin], analyticEvalD[ssE, rMin]};
     dSR = {analyticEval[tot, rMin], analyticEvalD[tot, rMin]} - dSS;
     rSS = threeRegionSolve[l, numericFunction[ssPsrc[Q, n, c, l, "real"], zVar], contourFunction[ssPsrc[Q, n, c, l, "continued"], zVar], dSS];
     srF = srReal[Q, n, c, l, lps]; rrF = rrReal[Q, n, c, l, lps];
     rSR = threeRegionSolve[l, srF[Automatic], srContour[Q, n, c, l, resSols], dSR];
     lSR = lineSolve[l, srF["in"], srF["out"], dSR];
     lRR = lineSolve[l, rrF["in"], rrF["out"], {0, 0}];
     Etot = (exteriorValue[rSS, #] + exteriorValue[rSR, #] + lRR["f"][#] - analyticEval[tot, #]) & /@ extPts;
     Ec = (exteriorValue[rSS, #] + lSR["f"][#] + lRR["f"][#] - analyticEval[tot, #]) & /@ extPts;
     dB = {Re[rSR["valR"]] - lSR["f"][aR], Re[rSR["derR"]] - lSR["d"][aR]};
     dec = decomp[Q, n, c, l, lps, nodes];
     out[l] = <|"scale" -> scale, "Etot" -> Max[Abs[Etot]]/scale, "Ec" -> Max[Abs[Ec]]/scale, "dB" -> dB, "a" -> dec["a"], "b" -> dec["b"],
       "aErr" -> homErr[l, dec["a"], extR, scale], "bErr" -> homErr[l, dec["b"], extR, scale], "dBerr" -> homErr[l, dB, extR, scale],
       "ssErr" -> Max[Abs[(exteriorValue[rSS, #] - analyticEval[ssE, #]) & /@ extPts]]/scale,
       "consistency" -> Max[Abs[dB - (dec["a"] + dec["b"])]], "consTol" -> 10^-6 Max[Abs[dB]] + 10^-16 scale|>],
    {l, lList}];
   out];

(* ================= main ================= *)
summary = {}; lpRows = {}; ampRows = {}; quadConv = {};
Do[
  SetCSParameters[casePars[Q]];
  Print["\n=== ", Q, " ==="];
  (* SS_P closed form: ODE check *)
  Do[Module[{e = ssPsol[Q, n, r0^-(n + 1), l]},
     csCheck[Q <> " n=" <> ToString[n] <> " l=" <> ToString[l] <> ": L_l(SS_P exact) = SS_P source (6 points, both sides)",
      Max[Table[Abs[evR[radialL[l, pick[e, If[r < r0, "in", "out"]], zVar], r] - evR[ssPsrc[Q, n, r0^-(n + 1), l, "real"], r]], {r, {1, 3, 9/2, 11/2, 7, 10}}]], 10^-25]],
   {n, nList}, {l, {0, 3, 6}}];
  (* main scan over (n, c) *)
  Do[
   Module[{n = st[[1]], c = st[[2]], res},
    res = runSetting[Q, n, c, lpMain, Range[0, lmax], "checkResidual" -> True];
    Do[Module[{x = res[l]},
       csCheck[Q <> " n=" <> ToString[n] <> " c=" <> ToString[c, InputForm] <> " l=" <> ToString[l] <> ": (a) real-part mismatch ~ 0 (scaled exterior response)", x["aErr"], 10^-16];
       csCheck[Q <> " n=" <> ToString[n] <> " c=" <> ToString[c, InputForm] <> " l=" <> ToString[l] <> ": contour - line SR at aR = (a) + (b)", x["consistency"], x["consTol"]];
       csCheck[Q <> " n=" <> ToString[n] <> " c=" <> ToString[c, InputForm] <> " l=" <> ToString[l] <> ": SS_P exterior error (scaled)", x["ssErr"], 10^-10];
       AppendTo[summary, {Q, n, c, l, x["Etot"], x["Ec"], x["aErr"], x["bErr"], x["dBerr"]}]], {l, 0, lmax}];
    Print["    n = ", n, ", c = ", ToString[c, InputForm], ":  max_l (b) = ", sci[Max[Table[res[l]["bErr"], {l, 0, lmax}]]], ",  max_l total = ",
     sci[Max[Table[res[l]["Etot"], {l, 0, lmax}]]], ",  max_l (c) = ", sci[Max[Table[res[l]["Ec"], {l, 0, lmax}]]], "   [", tick[], " s]"]],
   {st, settings}];
  (* quadrature convergence (one setting) *)
  Module[{d60 = decomp[Q, 3, r0^-4, 2, Range[0, lpMain], nodesFor[60]], d120 = decomp[Q, 3, r0^-4, 2, Range[0, lpMain], nodesFor[120]]},
   AppendTo[quadConv, Max[Abs[Join[d60["b"] - d120["b"], d60["a"] - d120["a"]]]]/Max[Abs[d60["b"]]]]];
  (* l'max convergence at n = 3 *)
  Do[Module[{res = runSetting[Q, 3, r0^-4, lpm, Range[0, lmax]]},
     Do[AppendTo[lpRows, {Q, lpm, l, res[l]["Etot"], res[l]["Ec"], res[l]["bErr"]}];
      csCheck[Q <> " l'max=" <> ToString[lpm] <> " l=" <> ToString[l] <> ": contour - line SR at aR = (a) + (b)", res[l]["consistency"], res[l]["consTol"]], {l, 0, lmax}];
     Print["    l'max = ", lpm, ": max_l (c) = ", sci[Max[Table[res[l]["Ec"], {l, 0, lmax}]]], ", max_l (b) = ", sci[Max[Table[res[l]["bErr"], {l, 0, lmax}]]], "   [", tick[], " s]"]],
   {lpm, lpScan}];
  (* amplitude independence at n = 3 *)
  Module[{base = runSetting[Q, 3, r0^-4, lpMain, Range[0, lmax]]},
   Do[SetCSParameters[Join[casePars[Q], <|"amp" -> amp|>]];
    Module[{res = runSetting[Q, 3, r0^-4, lpMain, Range[0, lmax]]},
     Do[csCheck[Q <> " amp=" <> ToString[amp, InputForm] <> " l=" <> ToString[l] <> ": contour SR shadow error dB(aR) unchanged", Max[Abs[res[l]["dB"] - base[l]["dB"]]], base[l]["consTol"]];
      AppendTo[ampRows, {Q, amp, l, Max[Abs[base[l]["dB"]]], Max[Abs[res[l]["dB"] - base[l]["dB"]]]}], {l, 0, lmax}]];
    Print["    amp = ", amp, " done [", tick[], " s]"], {amp, amps}];
   SetCSParameters[casePars[Q]]],
  {Q, {"Product", "GradDot"}}];
csCheck["quadrature converged (60 vs 120 nodes per side, relative, both cases)", Max[quadConv], 10^-10];

(* ================= (b) beyond l'max = 12: decomposition only (verified against the three-region solves up to l'max = 12) ================= *)
(* per-l' shadow contributions and their cumulative sum up to l' = lpProbe, n = 3, standard c; tube radius R = 1/4 and 1/2 *)
lpProbe = If[quick, 8, 30]; probeRows = {};
Do[SetCSParameters[casePars[Q]];
  Do[Block[{aL = r0 - R, aR = r0 + R, extR = Subdivide[r0 + R, rMax, 45]},
    Do[Module[{scale = Max[Abs[analyticEval[totalSol[Q, l], #] & /@ extPts]], per},
       per = Table[decomp[Q, 3, r0^-4, l, {lp}, nodesFor[40]]["b"], {lp, 0, lpProbe}];
       AppendTo[probeRows, {Q, R, l, Table[homErr[l, per[[lp + 1]], extR, scale], {lp, 0, lpProbe}],
         Table[homErr[l, Total[per[[1 ;; lpm + 1]]], extR, scale], {lpm, 0, lpProbe}]}]], {l, 0, lmax}]], {R, {1/4, 1/2}}];
  Print["    ", Q, ": l' probe done [", tick[], " s]"], {Q, {"Product", "GradDot"}}];
probeMax[Q_, R_, k_, idx_] := Max[Cases[probeRows, {Q, R, _, per_, cum_} :> {per, cum}[[k, idx + 1]]]];
probeLps = Select[{0, 4, 8, 12, 16, 20, 25, 30}, # <= lpProbe &];
md3 = StringRiffle[Join[{"| case | R | l' | (b) contribution of mode l' | cumulative (b), l'max = l' |", "|---|---|---|---|---|"},
    Flatten@Table["| " <> Q <> " | " <> ToString[R, InputForm] <> " | " <> ToString[lp] <> " | " <> sci[probeMax[Q, R, 1, lp]] <> " | " <> sci[probeMax[Q, R, 2, lp]] <> " |",
      {Q, {"Product", "GradDot"}}, {R, {1/2, 1/4}}, {lp, probeLps}]], "\n"];
growth = Table[{Q, R, N[(probeMax[Q, R, 1, lpProbe]/probeMax[Q, R, 1, lpProbe - 5])^(1/5), 4]}, {Q, {"Product", "GradDot"}}, {R, {1/2, 1/4}}];
Print["\n", md3, "\n\nper-l' growth factor of the (b) contribution near l' = ", lpProbe, " (max over l): ", growth];
Export[FileNameJoin[{resDir, "stage6_lprime_probe" <> tag <> ".csv"}],
  Prepend[Flatten[Table[{#[[1]], ToString[#[[2]], InputForm], #[[3]], lp, sci[#[[4, lp + 1]]], sci[#[[5, lp + 1]]]}, {lp, 0, lpProbe}] & /@ probeRows, 1],
   {"case", "R", "l", "l'", "b_contribution_of_mode_l'", "b_cumulative_lpmax=l'"}], "CSV"];

(* ================= tables and plot ================= *)
mx[Q_, n_, c_, k_] := Max[Cases[summary, {Q, n, c, __}][[All, k]]];
fitRate[pts_] := Module[{p = Select[pts, #[[2]] > 10^-17 &]}, If[Length[p] < 2, Missing[], 10^Coefficient[Fit[{#[[1]], Log10[#[[2]]]} & /@ p, {1, x}, x], x]]];
md1 = StringRiffle[Join[{"| case | n | c | (b) shadow | (a) mismatch | total error | (c) l' truncation (l'max = 12) | (b)/total |", "|---|---|---|---|---|---|---|---|"},
    Flatten@Table[Table["| " <> Q <> " | " <> ToString[st[[1]]] <> " | " <> ToString[st[[2]], InputForm] <> " | " <>
        StringRiffle[sci /@ {mx[Q, st[[1]], st[[2]], 8], mx[Q, st[[1]], st[[2]], 7], mx[Q, st[[1]], st[[2]], 5], mx[Q, st[[1]], st[[2]], 6], mx[Q, st[[1]], st[[2]], 8]/mx[Q, st[[1]], st[[2]], 5]}, " | "] <> " |",
       {st, settings}], {Q, {"Product", "GradDot"}}]], "\n"];
rates = Table[{Q, fitRate[Table[{n, mx[Q, n, r0^-(n + 1), 8]}, {n, nList}]]}, {Q, {"Product", "GradDot"}}];
linC = Table[{Q, Max[Table[Cases[summary, {Q, 3, 10 r0^-4, l, __}][[1, 8]]/Cases[summary, {Q, 3, r0^-4, l, __}][[1, 8]], {l, 0, lmax}]],
    Min[Table[Cases[summary, {Q, 3, 10 r0^-4, l, __}][[1, 8]]/Cases[summary, {Q, 3, r0^-4, l, __}][[1, 8]], {l, 0, lmax}]]}, {Q, {"Product", "GradDot"}}];
md2 = StringRiffle[Join[{"| case | l'max | max_l total | max_l (c) | max_l (b) |", "|---|---|---|---|---|"},
    Flatten@Table[Table["| " <> Q <> " | " <> ToString[lpm] <> " | " <> StringRiffle[sci /@ If[lpm == lpMain,
          {mx[Q, 3, r0^-4, 5], mx[Q, 3, r0^-4, 6], mx[Q, 3, r0^-4, 8]},
          {Max[Cases[lpRows, {Q, lpm, _, e_, _, _} :> e]], Max[Cases[lpRows, {Q, lpm, _, _, e_, _} :> e]], Max[Cases[lpRows, {Q, lpm, _, _, _, e_} :> e]]}], " | "] <> " |",
       {lpm, Append[lpScan, lpMain]}], {Q, {"Product", "GradDot"}}]], "\n"];
Print["\n", md1, "\n\nfitted decay of (b) per order of n (n = 1,3,5,7): ", rates, "\n(b)(10c)/(b)(c) at n = 3, max/min over l: ", linC, "\n\n", md2];
Export[FileNameJoin[{resDir, "stage6_summary" <> tag <> ".md"}],
  "# Stage 6: consistent 3D puncture\n\nScaled exterior errors (max over l = 0..6 of max|error|/max|analytic total| over [1,4.5] U [5.5,10]); Cauchy data at r = 1; K continued by Theta->H.\n\n" <> md1 <>
   "\n\nFitted decay factor of (b) per unit n: " <> ToString[rates] <> "\n\n(b)(10c)/(b)(c) at n = 3 (max, min over l): " <> ToString[linC] <> "\n\n## l'max convergence (n = 3)\n\n" <> md2 <> "\n\n## (b) per mode l' and cumulative, beyond l'max = 12 (decomposition only; n = 3)\n\n" <> md3 <>
   "\n\nPer-l' growth factor near l' = " <> ToString[lpProbe] <> ": " <> ToString[growth] <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage6_scan" <> tag <> ".csv"}], Prepend[{#[[1]], #[[2]], ToString[#[[3]], InputForm], #[[4]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]], sci[#[[9]]]} & /@ summary,
   {"case", "n", "c", "l", "total_err", "c_lprime_trunc_err", "a_mismatch", "b_shadow", "contour_minus_line_SR"}], "CSV"];
Export[FileNameJoin[{resDir, "stage6_lprime" <> tag <> ".csv"}], Prepend[{#[[1]], #[[2]], #[[3]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]]} & /@ lpRows, {"case", "l'max", "l", "total_err", "c_err", "b_err"}], "CSV"];
Export[FileNameJoin[{resDir, "stage6_amplitude" <> tag <> ".csv"}], Prepend[{#[[1]], ToString[#[[2]], InputForm], #[[3]], sci[#[[4]]], sci[#[[5]]]} & /@ ampRows, {"case", "amp", "l", "|dB(aR)| amp=1", "|dB(amp) - dB(1)|"}], "CSV"];
plt = ListLogPlot[Join[Table[Table[{n, mx[Q, n, r0^-(n + 1), 8]}, {n, nList}], {Q, {"Product", "GradDot"}}], Table[Table[{n, mx[Q, n, r0^-(n + 1), 6]}, {n, nList}], {Q, {"Product", "GradDot"}}],
    {Table[{n, 0.1^n}, {n, nList}]}], Joined -> True, PlotMarkers -> Automatic, Frame -> True, FrameLabel -> {"n", "scaled exterior error"},
   PlotStyle -> {Blue, Red, {Blue, Dashed}, {Red, Dashed}, {Gray, Dotted}},
   PlotLegends -> {"(b) shadow, 1/\[Rho]\.b2", "(b) shadow, 1/\[Rho]^4", "(c) \[ScriptL]' truncation, 1/\[Rho]\.b2", "(c) \[ScriptL]' truncation, 1/\[Rho]^4", "(R/\!\(\*SubscriptBox[\(r\), \(0\)]\)\!\(\*SuperscriptBox[\()\), \(n\)]\)"},
   BaseStyle -> {FontFamily -> "Times", FontSize -> 10}, ImageSize -> 425];
Export[FileNameJoin[{resDir, "stage6_error_vs_n" <> tag <> ".pdf"}], plt];
Print["\nStage 6 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage6_checks" <> tag <> ".csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
