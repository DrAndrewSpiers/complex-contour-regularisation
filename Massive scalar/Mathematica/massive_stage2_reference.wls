#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Massive scalar, Stage 2: the exact second-order reference.
   (Lap - m^2) Phi2 = (1/2) d_Z^2 (Phi1^2)  =>  Phi2 = (1/2) d_Z^2 chi with (Lap - m^2) chi = Phi1^2 = q^2 e^{-2 m rho}/rho^2 + 2 q e^{-m rho} PhiR/rho + PhiR^2.
   chi is not closed form in 3D (SR, RR), so the reference is built mode by mode:
     (1) closed-form modes of the three pieces of Phi1^2 (SS by the rho-substitution with exponential integrals; SR by Bessel couplings; RR by
         Gaunt couplings), checked against direct angular quadrature;
     (2) chi_j = Int G_j(r, s) [Phi1^2]_j(s) ds with the physical Green's function (NIntegrate, double-exponential rule, breakpoint at r0, 40 digits),
         chi_j'' and chi_j''' from the ODE, then the d_Z^2 ladder and the homogeneous correction A u_l + B v_l that meets the Robin conditions.
   Cross-checks: quadrature-rule independence; the SS piece against the 3D closed form chi_SS(rho) (exponential integrals) decomposed
   numerically; the SR and RR pieces against NDSolve real-line solves with the Robin conditions; m = 0 against the massless closed form chi.
   Writes results/massive_reference_<m>.mx (association used by Stages 3, 4) and CSV tables.
   Run:  wolframscript -file massive_stage2_reference.wls [smoke]                                                                      *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourMassive.wl"}]];
resDir = FileNameJoin[{ParentDirectory[dir], "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
smoke = MemberQ[$ScriptCommandLine, "smoke"];
csResetChecks[]; t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
Off[N::meprec]; Off[General::munfl]; Off[NIntegrate::ncvb];
mList = {1/5, 1, 0}; lmax = If[smoke, 2, 6];
r0 = 5; q = 1; c0 = 1; c1 = 1/10; rMin = 1; rMax = 10;
gIn = Subdivide[1, 7/2, 12]; gOut = Subdivide[13/2, 10, 14]; gCom = Join[gIn, gOut];
mStr[m_] := ToString[m, InputForm]; mTag[m_] := StringReplace[mStr[m], "/" -> "_"];
Print["=== Massive Stage 2: exact reference; m = ", mList, ", l <= ", lmax, " ==="];
refs = <||>; deltaRows = {}; srcRows = {}; xRows = {};

(* modes of a 3D function f[r, mu] at real r by angular quadrature *)
modeNum[f_, j_, r_] := Sqrt[Pi (2 j + 1)] NIntegrate[f[r, mu] LegendreP[j, mu], {mu, -1, 1}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 20];
rhoOf[r_, mu_] := Sqrt[r^2 + r0^2 - 2 r r0 mu];

Do[SetCMParameters[<|"m" -> m|>];
 Print["\n--- m = ", mStr[m], " ---"];
 (* (1) closed-form source modes against angular quadrature, both sides of r0 *)
 Module[{fSS, fSR, fRR},
  fSS[r_, mu_] := q^2 Exp[-2 m rhoOf[r, mu]]/rhoOf[r, mu]^2;
  fSR[r_, mu_] := 2 q Exp[-m rhoOf[r, mu]]/rhoOf[r, mu] phiR3D[r, mu];
  fRR[r_, mu_] := phiR3D[r, mu]^2;
  Do[Module[{br = If[r < r0, 1, 2], eSS, eSR, eRR},
     eSS = Abs[evR[logReal[ssUndiffPair["exact", j][[br]]], r] - modeNum[fSS, j, r]];
     eSR = Abs[evR[srUndiffPair["exact", j][[br]], r] - modeNum[fSR, j, r]];
     eRR = Abs[evR[rrUndiffPair["exact", j][[br]], r] - modeNum[fRR, j, r]];
     AppendTo[srcRows, {m, j, r, eSS, eSR, eRR}];
     csCheck["m=" <> mStr[m] <> " j=" <> ToString[j] <> " r=" <> ToString[r] <> ": closed-form modes of q^2 e^{-2m rho}/rho^2, 2q e^{-m rho} PhiR/rho, PhiR^2 = angular quadrature", eSS + eSR + eRR, 10^-22]],
   {j, If[smoke, {0, 2}, {0, 1, 3, 6, 8}]}, {r, {3, 49/10, 51/10, 8}}]];
 (* delta coefficients of the SS and SR source modes *)
 Do[Module[{dSS, dSR},
    dSS = (1/2) dzzDeltaCoeff[Function[j, ssUndiffPair["exact", j]], l];
    dSR = (1/2) dzzDeltaCoeff[Function[j, srUndiffPair["exact", j]], l];
    AppendTo[deltaRows, {m, l, dSS, dSR}];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": SS source delta coefficient = 4 Pi m q^2 Y_l(0)/r0^2 (from the -m q term of the singular field)", Abs[dSS - 4 Pi m q^2 Y0[l]/r0^2], 10^-25];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": SR source delta coefficient c_l = -4 Pi q PhiR(x0) Y_l(0)/r0^2 (m-independent)", Abs[dSR + 4 Pi q phiRAtParticle[] Y0[l]/r0^2], 10^-25]], {l, 0, lmax}];
 Print["    source modes and delta coefficients checked [", tick[], " s]"];
 (* (2) the reference, piece by piece *)
 Module[{ref = <||>},
  Do[ref[p] = exactReference[lmax, gCom, {p}]; Print["    reference piece ", p, " done [", tick[], " s]"], {p, {"SS", "SR", "RR"}}];
  ref["total"] = Association@Table[l -> <|"grid" -> gCom, "val" -> Sum[ref[p][l]["val"], {p, {"SS", "SR", "RR"}}], "der" -> Sum[ref[p][l]["der"], {p, {"SS", "SR", "RR"}}]|>, {l, 0, lmax}];
  refs[m] = ref;
  Do[csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": reference meets the Robin conditions at r = 1 and 10",
    Abs[robinResidual[l, ref["total"][l]["val"][[1]], ref["total"][l]["der"][[1]], 1]] + Abs[robinResidual[l, ref["total"][l]["val"][[-1]], ref["total"][l]["der"][[-1]], 10]], 10^-25], {l, 0, lmax}]];
 (* cross-check A: quadrature-rule independence (Gauss-Kronrod, j = 2, SS piece) *)
 Module[{alt, std, j = 2, nintSave},
  std = chiGreen[j, ssUndiffPair["exact", j], gCom];
  Block[{ContourMassive`Private`nint = Function[{f, a, b}, If[b <= a, 0, Module[{s}, NIntegrate[f[s], Evaluate[{s, Sequence @@ Union[{a, b}, Select[{r0}, a < # < b &]]}],
         Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule", "SingularityHandler" -> "IMT"}, WorkingPrecision -> 40, PrecisionGoal -> 28, MaxRecursion -> 40]]]]},
   alt = chiGreen[j, ssUndiffPair["exact", j], gCom]];
  AppendTo[xRows, {m, "A: DE vs Gauss-Kronrod/IMT, chi_2^SS and chi_2'", Max[Abs[std[[All, 1 ;; 2]] - alt[[All, 1 ;; 2]]]]/Max[Abs[std[[All, 1]]]]}];
  csCheck["m=" <> mStr[m] <> ": chi_2^SS by two NIntegrate rules agree (relative)", Max[Abs[std[[All, 1 ;; 2]] - alt[[All, 1 ;; 2]]]]/Max[Abs[std[[All, 1]]]], 10^-22]];
 (* cross-check B: SS piece against the 3D closed form chi_SS(rho) (m > 0 only: the closed form is written with exponential integrals) *)
 If[m != 0, Module[{rho, u, chiS, lhs, ser, d1, d2, src, vals, AB, err, mx, refSS = refs[m]["SS"]},
  chiS = q^2 ((Exp[m rho] ExpIntegralEi[-3 m rho] - Exp[-m rho] ExpIntegralEi[-m rho])/(2 m rho) - Log[3]/(2 m) Cosh[m rho]/rho);
  lhs = Simplify[D[chiS, {rho, 2}] + 2/rho D[chiS, rho] - m^2 chiS - q^2 Exp[-2 m rho]/rho^2];
  csCheck["m=" <> mStr[m] <> ": chi_SS(rho) = q^2 [(e^{m rho} Ei(-3 m rho) - e^{-m rho} Ei(-m rho))/(2 m rho) - ln3 cosh(m rho)/(2 m rho)] solves (Lap - m^2) chi = q^2 e^{-2 m rho}/rho^2 (symbolic)", Abs[N[lhs /. rho -> 7/3, 30]] + Abs[N[lhs /. rho -> 1/4, 30]], 10^-25];
  ser = Series[chiS, {rho, 0, 1}];
  csCheck["m=" <> mStr[m] <> ": chi_SS = q^2 ln rho + O(1): no 1/rho term (no point source at x0)", Abs[N[SeriesCoefficient[ser, -1], 30]] + Abs[N[Coefficient[Normal[ser], Log[rho]] - q^2, 30]], 10^-25];
  d1 = D[chiS, rho]; d2 = D[chiS, {rho, 2}];
  (* (1/2) d_Z^2 chi_SS = (1/2) [chi'' d^2/rho^2 + chi' (1/rho - d^2/rho^3)], d = r mu - r0 *)
  src = Function[{r, mu}, Evaluate[With[{dd = r mu - r0}, (1/2) ((d2 /. rho -> rhoOf[r, mu]) dd^2/rhoOf[r, mu]^2 + (d1 /. rho -> rhoOf[r, mu]) (1/rhoOf[r, mu] - dd^2/rhoOf[r, mu]^3))]]];
  Do[Module[{v, dv, ab},
     v = Table[modeNum[src, l, r], {r, gCom}];
     dv = Table[Sqrt[Pi (2 l + 1)] NIntegrate[D[src[rr, mu], rr] LegendreP[l, mu] /. rr -> r, {mu, -1, 1}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 20], {r, {1, 10}}];
     ab = physFix[l, v[[1]], dv[[1]], v[[-1]], dv[[2]]];
     err = Max[Abs[v + ab[[1]] (uFun[l, #] & /@ gCom) + ab[[2]] (vFun[l, #] & /@ gCom) - refSS[l]["val"]]]/Max[Abs[refSS[l]["val"]]];
     AppendTo[xRows, {m, "B: SS piece vs 3D closed form chi_SS(rho), l = " <> ToString[l], err}];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": SS reference (Green's quadrature of the mode source) = modes of (1/2) d_Z^2 chi_SS(rho) (scaled)", err, 10^-20]], {l, If[smoke, {0, 2}, {0, 1, 2, 4, 6}]}];
  Print["    cross-check B done [", tick[], " s]"]]];
 (* cross-check C: SR and RR pieces against NDSolve real-line solves (Robin conditions) *)
 Do[Module[{chis, vals, ab, err, pc = If[p === "SR", srUndiffPair["exact", #] &, rrUndiffPair["exact", #] &], refP = refs[m][p]},
    chis = Table[Module[{pr = pc[j], fIn, fOut, sol, S},
       fIn = numericFunction[pr[[1]], zVar]; fOut = numericFunction[pr[[2]], zVar];
       sol = realLinePhys[j, fIn, fOut];
       S[r_] := If[r < r0, fIn[r], fOut[r]];
       Table[With[{c = sol["f"][r], dc = sol["d"][r]}, Module[{d2, d3, dS},
          dS = If[r < r0, D[pr[[1]], zVar], D[pr[[2]], zVar]] /. zVar -> r;
          d2 = S[r] - 2 dc/r + (j (j + 1)/r^2 + m^2) c; d3 = dS - 2 d2/r + 2 dc/r^2 + (j (j + 1)/r^2 + m^2) dc - 2 j (j + 1) c/r^3; {c, dc, d2, d3}]], {r, gCom}]], {j, 0, lmax + 2}];
    Do[vals = Table[(1/2) dzzArr[chis[[All, k]], gCom[[k]]][[l + 1]], {k, Length[gCom]}];
     ab = physFix[l, vals[[1, 1]], vals[[1, 2]], vals[[-1, 1]], vals[[-1, 2]]];
     err = scaledError[vals[[All, 1]] + ab[[1]] (uFun[l, #] & /@ gCom) + ab[[2]] (vFun[l, #] & /@ gCom), refP[l]["val"]];       (* absolute if the piece vanishes (RR, l >= 5) *)
     AppendTo[xRows, {m, "C: " <> p <> " piece vs NDSolve (Robin), l = " <> ToString[l], err}];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": " <> p <> " reference = NDSolve real-line solve with the Robin conditions (scaled)", err, 10^-14], {l, 0, lmax}];
    Print["    cross-check C (", p, ") done [", tick[], " s]"]], {p, {"SR", "RR"}}];
 (* export *)
 Module[{ref = refs[m], tab},
  Export[FileNameJoin[{resDir, "massive_reference_m" <> mTag[m] <> ".mx"}], ref];
  tab = Flatten[Table[{mStr[m], l, p, N[gCom[[k]], 10], ToString[N[ref[p][l]["val"][[k]], 30], InputForm], ToString[N[ref[p][l]["der"][[k]], 30], InputForm]}, {p, {"SS", "SR", "RR", "total"}}, {l, 0, lmax}, {k, Length[gCom]}], 2];
  Export[FileNameJoin[{resDir, "massive_reference_m" <> mTag[m] <> ".csv"}], Prepend[tab, {"m", "l", "piece", "r", "Phi2_l", "dPhi2_l"}], "CSV"]],
 {m, mList}];

(* ================= m = 0: the massless closed form ================= *)
Print["\n--- m = 0 (massless closed form) ---"];
SetCMParameters[<|"m" -> 0|>];
Module[{chiML, ref0, rr, mu, Z, rho, chiF, dchiF, src0, errs = {}},
  ref0 = refs[0];
  (* PRD App. A: chi = q^2 ln rho + q (c0 + c1 r0) rho + (1/2) q c1 (Z - r0) rho + c0^2 r^2/6 + c0 c1 r^2 Z/5 + c1^2 Z^4/12 *)
  chiML = q^2 Log[rho] + q (c0 + c1 r0) rho + (1/2) q c1 (Z - r0) rho + c0^2 rr^2/6 + c0 c1 rr^2 Z/5 + c1^2 Z^4/12 /. {Z -> rr mu, rho -> Sqrt[rr^2 + r0^2 - 2 rr r0 mu]};
  chiF = Function[{r, m2}, Evaluate[chiML /. {rr -> r, mu -> m2}]]; dchiF = Function[{r, m2}, Evaluate[D[chiML, rr] /. {rr -> r, mu -> m2}]];
  Module[{chis, vals, ab, err},
     chis = Table[Module[{pr = ssUndiffPair["exact", j] + srUndiffPair["exact", j] + rrUndiffPair["exact", j], S, dS},
        S[r_] := evR[logReal[pr[[If[r < r0, 1, 2]]]], r]; dS[r_] := evR[D[logReal[pr[[If[r < r0, 1, 2]]]], zVar], r];
        Table[With[{c = modeNum[chiF, j, r], dc = modeNum[dchiF, j, r]}, Module[{d2, d3},
           d2 = S[r] - 2 dc/r + (j (j + 1)/r^2) c; d3 = dS[r] - 2 d2/r + 2 dc/r^2 + (j (j + 1)/r^2) dc - 2 j (j + 1) c/r^3; {c, dc, d2, d3}]], {r, gCom}]], {j, 0, lmax + 2}];
     Do[vals = Table[(1/2) dzzArr[chis[[All, k]], gCom[[k]]][[l + 1]], {k, Length[gCom]}];
      ab = physFix[l, vals[[1, 1]], vals[[1, 2]], vals[[-1, 1]], vals[[-1, 2]]];
      err = Max[Abs[vals[[All, 1]] + ab[[1]] (uFun[l, #] & /@ gCom) + ab[[2]] (vFun[l, #] & /@ gCom) - ref0["total"][l]["val"]]]/Max[Abs[ref0["total"][l]["val"]]];
      AppendTo[xRows, {0, "m = 0: Green's quadrature vs massless closed form chi, l = " <> ToString[l], err}];
      csCheck["m=0 l=" <> ToString[l] <> ": reference machinery at m = 0 = modes of the massless closed form (1/2) d_Z^2 chi (scaled)", err, 10^-18], {l, 0, lmax}]];
  Print["    m = 0 closed-form comparison done [", tick[], " s]"];
  (* no continuity check at small m > 0: the closed forms with Exp[-c rho] cancel like 1/c^(2j+p+2) and are not usable at 50 digits below
     m ~ 0.1 (an m = 1/1000 and an m = 1/20 attempt both hung in NIntegrate); the m -> 0 limit is checked exactly at m = 0 above *)
  SetCMParameters[<|"m" -> 0|>]];
Do[csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": SS + SR + RR reference pieces recorded",
  If[Length[refs[m]["total"][l]["val"]] == Length[gCom], 0, 1], 0], {m, {1/5, 1, 0}}, {l, {0}}];

(* ================= tables and figure ================= *)
mdD = StringRiffle[Join[{"| m | l | SS delta coefficient | SR delta coefficient c_l |", "|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[NumberForm[N[#[[3]], 8], 8]] <> " | " <> ToString[NumberForm[N[#[[4]], 8], 8]] <> " |") & /@ deltaRows], "\n"];
mdS = StringRiffle[Join[{"| m | j | r | SS | SR | RR |", "|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]], InputForm] <> " | " <> StringRiffle[sci /@ #[[4 ;; 6]], " | "] <> " |") & /@ srcRows], "\n"];
mdX = StringRiffle[Join[{"| m | cross-check | scaled difference |", "|---|---|---|"}, ("| " <> mStr[#[[1]]] <> " | " <> #[[2]] <> " | " <> sci[#[[3]]] <> " |") & /@ xRows], "\n"];
mdV = StringRiffle[Join[{"| m | l | max over exterior grid of abs(Phi2_l) | SS | SR | RR |", "|---|---|---|---|---|---|"},
    Flatten@Table["| " <> mStr[m] <> " | " <> ToString[l] <> " | " <> StringRiffle[sci[Max[Abs[refs[m][#][l]["val"]]]] & /@ {"total", "SS", "SR", "RR"}, " | "] <> " |", {m, {1/5, 1, 0}}, {l, 0, lmax}]], "\n"];
Print["\n", mdD, "\n\n", mdX, "\n\n", mdV];
Export[FileNameJoin[{resDir, "massive_stage2_summary.md"}], "# Massive Stage 2: exact reference\n\n## delta coefficients of the source modes\n\n" <> mdD <> "\n\n## closed-form source modes against angular quadrature (absolute differences)\n\n" <> mdS <>
   "\n\n## cross-checks of the reference\n\n" <> mdX <> "\n\n## size of the reference modes (scale of the errors quoted later)\n\n" <> mdV <> "\n", "Text"];
Export[FileNameJoin[{resDir, "massive_stage2_delta_coefficients.csv"}], Prepend[{mStr[#[[1]]], #[[2]], N[#[[3]], 15], N[#[[4]], 15]} & /@ deltaRows, {"m", "l", "SS_delta", "SR_delta"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage2_crosschecks.csv"}], Prepend[{mStr[#[[1]]], #[[2]], sci[#[[3]]]} & /@ xRows, {"m", "check", "scaled_difference"}], "CSV"];
(* figure: the reference modes l = 0, 2, 4 for m = 0, 1/5, 1 (CSV + PDF) *)
$fs = {FontFamily -> "Times", FontSize -> 9};
Module[{cols = <|0 -> Black, 1/5 -> RGBColor[0.1, 0.35, 0.75], 1 -> RGBColor[0.85, 0.45, 0.1]|>, pl, data},
  data = Flatten[Table[{mStr[m], l, N[gCom[[k]], 6], N[refs[m]["total"][l]["val"][[k]], 8]}, {m, {0, 1/5, 1}}, {l, {0, 2, 4}}, {k, Length[gCom]}], 2];
  Export[FileNameJoin[{figDir, "fig_reference_modes.csv"}], Prepend[data, {"m", "l", "r", "Phi2_l"}], "CSV"];
  pl = ListPlot[Flatten[Table[{Table[{gCom[[k]], refs[m]["total"][l]["val"][[k]]}, {k, Length[gIn]}], Table[{gCom[[k]], refs[m]["total"][l]["val"][[k]]}, {k, Length[gIn] + 1, Length[gCom]}]}, {l, {0, 2, 4}}, {m, {0, 1/5, 1}}], 2],
    Joined -> True, PlotStyle -> Flatten[Table[{Directive[cols[m], Which[l == 0, Dashing[{}], l == 2, Dashed, True, Dotted], AbsoluteThickness[1.2]], Directive[cols[m], Which[l == 0, Dashing[{}], l == 2, Dashed, True, Dotted], AbsoluteThickness[1.2]]}, {l, {0, 2, 4}}, {m, {0, 1/5, 1}}], 2],
    Frame -> True, Axes -> False, FrameLabel -> {"r", "\[CapitalPhi]\!\(\*SubscriptBox[\(2\), \(\[ScriptL]\)]\)(r)"}, BaseStyle -> $fs, ImageSize -> 300, PlotRange -> All,
    Prolog -> {GrayLevel[0.9], Rectangle[{9/2, -10^3}, {11/2, 10^3}]},
    Epilog -> {Inset[Style["black m = 0, blue m = 1/5, orange m = 1; solid \[ScriptL] = 0, dashed 2, dotted 4", $fs], Scaled[{0.97, 0.95}], {Right, Top}]}];
  Export[FileNameJoin[{figDir, "fig_reference_modes.pdf"}], pl]];
Print["\nStage 2 wall time: ", tick[], " s"];
csSummary[];
exportChecks[FileNameJoin[{resDir, "massive_stage2_checks.csv"}]];
