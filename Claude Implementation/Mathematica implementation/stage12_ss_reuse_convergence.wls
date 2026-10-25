#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 12 (Dzz case, Q = (1/2) d_z^2 (A B)).
   (A) Reuse of the hybrid technique for SS (note Sec. 14.5): exact puncture, l = 0..6.
       SS_l = singular part [A/Dr^2 + B/Dr + C0 ln|Dr| (+ C1 Dr ln|Dr| + C2 Dr^2 ln|Dr|)], continued analytically,
            + remainder, known only as samples on a real grid (spacing h) and Cauchy-continued:
              s~(z) = (1/(Pi i)) Int_I s(t)/(t - z) dt,  I = [4, 6] > tube [4.5, 5.5].
       The remainder is represented by local degree-p Lagrange interpolants (one-sided at r0); the Cauchy integral of
       the interpolant is exact cell by cell near Re z and 16-point Gauss elsewhere.  Exterior SS from the Green's-function
       contour quadrature (equivalent to the contour solve by Lemma 1; one contour ODE solve cross-checks it).
   (B) Hybrid H1 (k = 0) convergence for the Stage 6 puncture, n = 1, 3, 5, 7, l <= 6: error against l'max, per-l'
       contributions, error against r, and tail acceleration (power-law and Upton's Lambda_{l,0,n}, arXiv:2508.00087 eq. 173).
   Run:  wolframscript -file stage12_ss_reuse_convergence.wls [A] [B]   (default: both)                                  *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
s9 = Import[FileNameJoin[{dir, "stage9_dzz.wls"}], "Text"];
defsA = StringTake[s9, {1, StringPosition[s9, "(* ================= checks of the d_z machinery"][[1, 1]] - 1}];
defsA = StringReplace[StringReplace[defsA, "dir = DirectoryName[$InputFileName];" -> "", 1], "Print[\"=== Stage 9" -> "Hold[\"x", 1];
defsB = StringTake[s9, {StringPosition[s9, "(* ================= TEST 2"][[1, 1]], StringPosition[s9, "rows9 = {};"][[1, 1]] - 1}];
defsB = StringReplace[defsB, "Print[\"\\n--- Test 2" -> "Hold[\"x", 1];
ToExpression[defsA]; ToExpression[defsB];                 (* Stage 6-9 definitions; SetCSParameters[casePars["Product"]] *)
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
args = Rest[$ScriptCommandLine]; phys = MemberQ[args, "phys"]; args = DeleteCases[args, "phys"];
{doA, doB, doC, doD} = (args === {} || MemberQ[args, #]) & /@ {"A", "B", "C", "D"};
If[phys, doB = False];          (* physical boundary conditions (psi' = l psi/r at 1, -(l+1) psi/r at 10) for A, C, D; B has its own script, stage12b *)
tagR = If[phys, "stage12phys_", "stage12_"];
physFix[sol_, l_] := Module[{B = (analyticEvalD[sol, 1] - l analyticEval[sol, 1])/(2 l + 1), A = -(analyticEvalD[sol, 10] + (l + 1) analyticEval[sol, 10]/10)/((2 l + 1) 10^(l - 1))},
   Function[r, analyticEval[sol, r] + A r^l + B r^(-l - 1)]];
If[(doC || doD) && ! doA, doA = True];               (* C and D reuse the A set-up *)
$fs = {FontFamily -> "Times", FontSize -> 9};
Print["=== Stage 12: SS reuse of the hybrid technique (A: ", doA, "), H1 convergence (B: ", doB, ") ==="];

(* ======================================================================================= *)
(* (A) SS: analytic singular part + Cauchy-continued numerical remainder                    *)
(* ======================================================================================= *)
If[doA, Module[{},
  Print["\n--- (A) SS reuse of the hybrid technique ---"];
  lA = Range[0, 6]; hList = {1/50, 1/100, 1/200}; epsList = {10^-12, 10^-8}; vNames = {"A", "AB", "ABC", "ABC+"};
  masks = <|"A" -> {1, 0, 0, 0, 0}, "AB" -> {1, 1, 0, 0, 0}, "ABC" -> {1, 1, 1, 0, 0}, "ABC+" -> {1, 1, 1, 1, 1}|>;
  sRealE[l_] := ssSourceModeExpr["Dzz", l, zVar, "real"] /. Log[x_^2] :> 2 Log[Abs[x]];
  (* local expansion at r0: S = A/Dr^2 + B/Dr + (C0 + C1 Dr + C2 Dr^2) ln|Dr| + analytic; from both sides *)
  sCoefBoth[l_] := Module[{e = sRealE[l], u, sp, sm, cp, cm},
     sp = Normal[Series[(e /. zVar -> r0 + u) /. Abs[u] -> u, {u, 0, 2}]];
     sm = Normal[Series[(e /. zVar -> r0 - u) /. Abs[-u] -> u /. Abs[u] -> u, {u, 0, 2}]];
     cp = {Coefficient[sp, u, -2], Coefficient[sp, u, -1], Coefficient[sp, Log[u]] /. u -> 0, Coefficient[Coefficient[sp, Log[u]], u, 1], Coefficient[Coefficient[sp, Log[u]], u, 2]};
     cm = {Coefficient[sm, u, -2], -Coefficient[sm, u, -1], Coefficient[sm, Log[u]] /. u -> 0, -Coefficient[Coefficient[sm, Log[u]], u, 1], Coefficient[Coefficient[sm, Log[u]], u, 2]};
     {cp, cm}];
  sing[cf_, mk_, d_, lg_] := (mk cf) . {1/d^2, 1/d, lg, d lg, d^2 lg};
  (* contour Gauss nodes in a (40 pieces x 30), z = y(a), weights y'(a) w *)
  yA[a_] := a + I (1 - ((a - r0)/(1/2))^2)^6;
  cNodesZ = Module[{br = Subdivide[r0 - 1/2, r0 + 1/2, 40], a}, Flatten[Table[({N[yA[#[[1]]], 40], N[D[yA[a], a] /. a -> #[[1]], 40], #[[2]]} & /@
         GaussianQuadratureWeights[30, br[[k]], br[[k + 1]], 40]), {k, Length[br] - 1}], 1]];
  zs = cNodesZ[[All, 1]]; zsM = N[zs]; dzwHP = cNodesZ[[All, 2]] cNodesZ[[All, 3]]; dzw = N[dzwHP];
  {ug, wg} = Transpose[N[GaussianQuadratureWeights[16, -1, 1, 30]]];
  vinvP[p_] := vinvP[p] = Table[N[Inverse[Table[If[k == 0, 1, (2 (i - d))^k], {i, 0, p}, {k, 0, p}]]], {d, 0, p}];
  (* local degree-p Lagrange interpolant per cell (one-sided at r0 and at the ends of I) *)
  cauchySetup[rem_, h_, p_] := Module[{nC = Length[rem], half = Length[rem]/2, tc, coef, upow = Table[ug^k, {k, 0, p}]},
     tc = N[4 + (Range[nC] - 1/2) h];
     coef = Table[Module[{lo = If[c <= half, 1, half + 1], hi = If[c <= half, half, nC], s0},
         s0 = Clip[c - Floor[p/2], {lo, hi - p}]; vinvP[p][[c - s0 + 1]] . rem[[s0 ;; s0 + p]]], {c, nC}];
     <|"tc" -> tc, "coef" -> coef, "h" -> h, "p" -> p, "upow" -> upow, "tG" -> Flatten[Outer[Plus, tc, (h/2) ug]], "fG" -> Flatten[coef . upow], "wG" -> Flatten[ConstantArray[(h/2) wg, nC]]|>];
  (* exact Int_{-1}^{1} P(u)/(u - zeta) du for a cell polynomial with coefficients a *)
  exactCell[a_, zeta_] := Module[{p = Length[a] - 1, pz = Fold[#1 zeta + #2 &, 0, Reverse[a]]},
     pz (Log[1 - zeta] - Log[-1 - zeta]) + Sum[a[[k + 1]] Sum[If[EvenQ[j], 2/(j + 1) zeta^(k - 1 - j), 0], {j, 0, k - 1}], {k, 1, p}]];
  (* cells with |zeta| <= 3 exactly, all others by 16-point Gauss; near cells are masked out of the Gauss sum (not subtracted:
     for z almost on the axis a Gauss term can be ~1/Im z and its subtraction would cancel catastrophically) *)
  cauchyAt2[cs0_, z0_] := Module[{z = If[Im[z0] == 0, z0 + I 10^-30, z0], zeta, nearQ, far, corr = 0},
     zeta = (z - cs0["tc"])/(cs0["h"]/2); nearQ = UnitStep[3 - Abs[zeta]];
     far = Total[(1 - Flatten[Outer[Times, nearQ, ConstantArray[1, Length[ug]]]]) cs0["wG"] cs0["fG"]/(cs0["tG"] - z)];
     Do[corr += exactCell[cs0["coef"][[c]], zeta[[c]]], {c, Flatten[Position[nearQ, 1]]}];
     (far + corr)/(Pi I)];
  (* contour-source difference (hybrid - closed form) at the nodes, and the exterior error it causes (Green's function, zero data at r = 1) *)
  dSource[l_, rem_, h_, mk_, p_] := Module[{cs = cauchySetup[rem, h, p]},
     N[sing[cfA[l], mk, zsM - r0, Log[zsM - r0]]] + (cauchyAt2[cs, #] & /@ zsM) - sCfM[l]];
  extErr[l_, ds_] := Module[{dM = {Re[Total[zsM^(1 - l) ds dzw]], Re[Total[zsM^(l + 2) ds dzw]]}},
     If[phys, Table[If[r < r0, -dM[[1]] r^l, -dM[[2]] r^(-l - 1)]/(2 l + 1), {r, gCom}],       (* free-space Green's function, both exterior regions *)
      (#^l dM[[1]] - #^(-l - 1) dM[[2]])/(2 l + 1) & /@ gOut]];
  aRows = {}; nRows = {}; cfRows = {};
  Do[Module[{e2 = sCoefBoth[l], fCf, fRe, solSS, scale, data, hom, MCf, psiCf, tail, ex},
     csCheck["(A) l=" <> ToString[l] <> ": singular coefficients A, B, C0, C1, C2 agree from both sides of r0", Max[Abs[N[e2[[1]] - e2[[2]], 30]]], 10^-25];
     cfA[l] = N[e2[[1]], 60];
     fCf = Function @@ {{zVar}, N[ssSourceModeExpr["Dzz", l, zVar, "continued"], 50]};
     fRe = Function @@ {{zVar}, N[sRealE[l], 50]};
     sCfHP = fCf /@ zs; sCfM[l] = N[sCfHP];
     (* closed-form contour against exact (high precision) *)
     solSS = pieceAnalyticSolution["Dzz", "SS", l]; scale = If[phys, Max[Abs[physFix[solSS, l] /@ gCom]], Max[Abs[analyticEval[solSS, #] & /@ gCom]]]; scaleA[l] = scale;
     data = {analyticEval[solSS, 1], analyticEvalD[solSS, 1]};
     hom = Module[{A, B}, (A #^l + B #^(-l - 1)) & /. First@Solve[{A + B == data[[1]], l A - (l + 1) B == data[[2]]}, {A, B}]];
     MCf = Table[Re[Total[zs^w sCfHP dzwHP]], {w, {1 - l, l + 2}}] +
       Table[Module[{s}, NIntegrate[s^w fRe[s], {s, 1, aL}, WorkingPrecision -> 40, PrecisionGoal -> 25] + NIntegrate[s^w fRe[s], {s, aR, 13/2}, WorkingPrecision -> 40, PrecisionGoal -> 25]], {w, {1 - l, l + 2}}];
     tail = Association@Table[r -> Table[Module[{s}, NIntegrate[s^w fRe[s], {s, 13/2, r}, WorkingPrecision -> 40, PrecisionGoal -> 25]], {w, {1 - l, l + 2}}], {r, gOut}];
     psiCf = (hom[#] + (#^l (MCf[[1]] + tail[#][[1]]) - #^(-l - 1) (MCf[[2]] + tail[#][[2]]))/(2 l + 1)) & /@ gOut;
     ex = analyticEval[solSS, #] & /@ gOut;
     If[phys, Module[{Mt = Table[Re[Total[zs^w sCfHP dzwHP]], {w, {1 - l, l + 2}}], ni},
        ni[w_, lo_, hi_] := If[hi <= lo, 0, Module[{s}, NIntegrate[s^w fRe[s], {s, lo, hi}, WorkingPrecision -> 40, PrecisionGoal -> 25]]];
        psiCf = Table[-(r^(-l - 1) (ni[l + 2, 1, Min[r, aL]] + If[r > aR, Mt[[2]] + ni[l + 2, aR, r], 0]) +
              r^l (If[r < aL, ni[1 - l, r, aL] + Mt[[1]] + ni[1 - l, aR, 10], ni[1 - l, r, 10]]))/(2 l + 1), {r, gCom}];
        ex = physFix[solSS, l] /@ gCom]];
     cfExact[l] = psiCf - ex;
     AppendTo[cfRows, {l, Max[Abs[cfExact[l]]]/scale}];
     csCheck["(A) l=" <> ToString[l] <> ": closed-form contour SS (Green's-function quadrature) = exact (scaled)", Max[Abs[cfExact[l]]]/scale, 10^-20];
     (* remainder samples and the hybrid variants *)
     Do[Module[{tcs = 4 + (Range[2/h] - 1/2) h, sV, xi},
        sV = fRe /@ tcs; SeedRandom[1000 l + Round[1/h]]; xi = RandomVariate[NormalDistribution[], Length[tcs]];
        Do[Module[{rem = N[sV - (sing[cfA[l], masks[v], # - r0, Log[Abs[# - r0]]] & /@ tcs)], e},
           Do[If[p == 7 || v === "ABC+",
             e = extErr[l, dSource[l, rem, h, masks[v], p]];
             AppendTo[aRows, {l, v, h, p, Max[Abs[e]]/scale, Max[Abs[e + cfExact[l]]]/scale}]], {p, {7, 3}}];
           If[MemberQ[{"ABC", "ABC+"}, v],
            Do[Module[{remN = If[model === "rel. remainder", rem (1 + eps xi), rem + eps xi Abs[N[sV]]], en},
               en = extErr[l, dSource[l, remN, h, masks[v], 7]];
               AppendTo[nRows, {l, v, h, model, eps, Max[Abs[en]]/scale, Max[Abs[en - extErr[l, dSource[l, rem, h, masks[v], 7]]]]/scale}]],
             {model, {"rel. remainder", "rel. full SS"}}, {eps, epsList}]]], {v, vNames}]], {h, hList}];
     Print["    l=", l, " [", tick[], " s]"]], {l, lA}];
  (* Lemma-1 cross-check of the Cauchy machinery (l = 2, ABC, h = 1/100): the contour quadrature of G(aR, z) dS(z) must equal the
     real-line integral of G(aR, x) Re dS(x + i0) over the tube, and Re dS(x + i0) = interpolant - exact remainder (pointwise, bounded) *)
  Module[{l = 2, h = 1/100, tcs, fRe2, rem, cs, quad, line, sel, xg, g0a, g1a, g0x, g1x, exRem, x},
    tcs = 4 + (Range[2/h] - 1/2) h;
    fRe2 = Function @@ {{zVar}, N[sRealE[l], 80]};
    rem = N[(fRe2 /@ tcs) - (sing[cfA[l], masks["ABC"], # - r0, Log[Abs[# - r0]]] & /@ tcs)];
    cs = cauchySetup[rem, h, 7];
    g0a = zsM^2/(2 l + 1) (aR^l zsM^(-l - 1) - zsM^l aR^(-l - 1)); g1a = zsM^2/(2 l + 1) (l aR^(l - 1) zsM^(-l - 1) + (l + 1) zsM^l aR^(-l - 2));
    quad = Re[{Total[g0a dSource[l, rem, h, masks["ABC"], 7] dzw], Total[g1a dSource[l, rem, h, masks["ABC"], 7] dzw]}];
    sel = Flatten[Position[cs["tG"], x_ /; aL < x < aR]]; xg = cs["tG"][[sel]];
    g0x = xg^2/(2 l + 1) (aR^l xg^(-l - 1) - xg^l aR^(-l - 1)); g1x = xg^2/(2 l + 1) (l aR^(l - 1) xg^(-l - 1) + (l + 1) xg^l aR^(-l - 2));
    (* interpolant part: Gauss per cell (exact for the polynomial); exact-remainder part: adaptive quadrature at 40 digits (Dr ln Dr at r0) *)
    line = {Total[cs["wG"][[sel]] g0x cs["fG"][[sel]]], Total[cs["wG"][[sel]] g1x cs["fG"][[sel]]]} -
      Table[N[NIntegrate[gk (fRe2[x] - sing[cfA[l], masks["ABC"], x - r0, Log[Abs[x - r0]]]), {x, aL, r0, aR}, WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 40]],
       {gk, {x^2/(2 l + 1) (aR^l x^(-l - 1) - x^l aR^(-l - 1)), x^2/(2 l + 1) (l aR^(l - 1) x^(-l - 1) + (l + 1) x^l aR^(-l - 2))}}];
    odeCheck = {line, quad};
    csCheck["(A) Lemma 1 for the hybrid source (l=2, ABC, h=1/100): contour quadrature = real-line integral of (interpolant - remainder) (relative)", Max[Abs[line - quad]]/Max[Abs[quad]], 10^-5]];
  (* checks and tables *)
  mxA[v_, h_, p_, k_] := Max[Cases[aRows, {_, v, h, p, a_, b_} :> {a, b}[[k]]]];
  csCheck["(A) {A,B,C0}, p = 7: error falls at least as h^2.5 from h = 1/100 to 1/200 (shortfall in the order)", Max[0, 5/2 - Log[2, mxA["ABC", 1/100, 7, 1]/mxA["ABC", 1/200, 7, 1]]], 0];
  csCheck["(A) {A,B,C0,C1,C2}, h = 1/200, p = 7: exterior SS = closed-form contour (max over l, scaled)", mxA["ABC+", 1/200, 7, 1], 10^-10];
  Do[csCheck["(A) negative control: subtracting only {" <> v <> "} is first order and over 100x worse than {A,B,C0} at h = 1/200", If[mxA[v, 1/200, 7, 1] > 100 mxA["ABC", 1/200, 7, 1] && Abs[Log[2, mxA[v, 1/100, 7, 1]/mxA[v, 1/200, 7, 1]] - 1] < 0.2, 0, 1], 0], {v, {"A", "AB"}}];
  mdA1 = StringRiffle[Join[{"| subtracted | p | h = 1/50 | 1/100 | 1/200 | order (1/100 -> 1/200) |", "|---|---|---|---|---|---|"},
      Flatten@Table[If[p == 7 || v === "ABC+", {"| {" <> v <> "} | " <> ToString[p] <> " | " <> StringRiffle[sci[mxA[v, #, p, 1]] & /@ hList, " | "] <> " | " <>
           ToString[NumberForm[Log[2, mxA[v, 1/100, p, 1]/mxA[v, 1/200, p, 1]], 3]] <> " |"}, {}], {v, vNames}, {p, {7, 3}}]], "\n"];
  mxN[v_, h_, model_, eps_, k_] := Max[Cases[nRows, {_, v, h, model, eps, a_, b_} :> {a, b}[[k]]]];
  mdA2 = StringRiffle[Join[{"| subtracted | h | noise model | eps | total error | noise part |", "|---|---|---|---|---|---|"},
      Flatten@Table["| {" <> v <> "} | " <> ToString[h, InputForm] <> " | " <> model <> " | " <> sci[N[eps]] <> " | " <> sci[mxN[v, h, model, eps, 1]] <> " | " <> sci[mxN[v, h, model, eps, 2]] <> " |",
        {v, {"ABC", "ABC+"}}, {h, hList}, {model, {"rel. remainder", "rel. full SS"}}, {eps, epsList}]], "\n"];
  mdA3 = "closed-form contour vs exact (max over l, scaled): " <> sci[Max[cfRows[[All, 2]]]] <> "; Lemma-1 check at aR (l = 2, ABC, h = 1/100), contour vs real line: " <>
    sci[Max[Abs[odeCheck[[1]] - odeCheck[[2]]]]/Max[Abs[odeCheck[[2]]]]] <> " relative.";
  Print["\n", mdA1, "\n\n", mdA2, "\n\n", mdA3];
  Export[FileNameJoin[{resDir, tagR <> "ss_reuse.csv"}], Prepend[{#[[1]], #[[2]], ToString[#[[3]], InputForm], #[[4]], sci[#[[5]]], sci[#[[6]]]} & /@ aRows, {"l", "subtracted", "h", "p", "err_vs_closed_form_contour", "err_vs_exact"}], "CSV"];
  Export[FileNameJoin[{resDir, tagR <> "ss_noise.csv"}], Prepend[{#[[1]], #[[2]], ToString[#[[3]], InputForm], #[[4]], sci[N[#[[5]]]], sci[#[[6]]], sci[#[[7]]]} & /@ nRows, {"l", "subtracted", "h", "noise_model", "eps", "err_vs_cf", "noise_part"}], "CSV"];
  Export[FileNameJoin[{resDir, tagR <> "A_summary.md"}], "# Stage 12 (A): SS reuse of the hybrid technique\n\nScaled exterior SS error (max over l = 0..6) against the closed-form contour result.\n\n" <> mdA1 <>
     "\n\n## noise on the remainder samples\n\n" <> mdA2 <> "\n\n" <> mdA3 <> "\n", "Text"];
  Print["    (A) done [", tick[], " s]"]]];

(* ======================================================================================= *)
(* (C) subtraction-free: radial factor F = (r - r0)^2 SS, Cauchy-continued, pinned          *)
(* ======================================================================================= *)
If[doC, Module[{},
  Print["\n--- (C) radial factor F = (r - r0)^2 SS ---"];
  (* continuous piecewise interpolant: pieces between consecutive samples (and to the ends of I), centred stencils of p+1 samples *)
  vinvGen[nodes_] := vinvGen[nodes] = N[Inverse[Table[If[k == 0, 1, u^k], {u, nodes}, {k, 0, Length[nodes] - 1}]]];
  pieceSetup[vals_, h_, p_] := Module[{nS = Length[vals], ts, brk, ms, ws, coef, upow = Table[ug^k, {k, 0, p}]},
     ts = 4 + (Range[nS] - 1/2) h; brk = Join[{4}, ts, {6}]; ms = (Most[brk] + Rest[brk])/2; ws = (Rest[brk] - Most[brk])/2;
     coef = Table[Module[{s0 = Clip[k - (p + 1)/2, {1, nS - p}]}, vinvGen[(ts[[s0 ;; s0 + p]] - ms[[k]])/ws[[k]]] . vals[[s0 ;; s0 + p]]], {k, Length[ms]}];
     <|"tc" -> N[ms], "w" -> N[ws], "coef" -> coef, "upow" -> upow, "tG" -> Flatten[MapThread[#1 + #2 ug &, {N[ms], N[ws]}]], "fG" -> Flatten[coef . upow], "wG" -> Flatten[(# wg) & /@ N[ws]]|>];
  derivSetup[cs_] := Module[{p = Length[cs["coef"][[1]]] - 1, cD},
     cD = MapThread[Append[Table[(k + 1) #1[[k + 2]], {k, 0, p - 1}]/#2, 0] &, {cs["coef"], cs["w"]}];
     Join[cs, <|"coef" -> cD, "fG" -> Flatten[cD . cs["upow"]]|>]];
  cauchyAtP[cs0_, z0_] := Module[{z = If[Im[z0] == 0, z0 + I 10^-30, z0], zeta, nearQ, far, corr = 0},
     zeta = (z - cs0["tc"])/cs0["w"]; nearQ = UnitStep[3 - Abs[zeta]];
     far = Total[(1 - Flatten[Outer[Times, nearQ, ConstantArray[1, Length[ug]]]]) cs0["wG"] cs0["fG"]/(cs0["tG"] - z)];
     Do[corr += exactCell[cs0["coef"][[c]], zeta[[c]]], {c, Flatten[Position[nearQ, 1]]}];
     (far + corr)/(Pi I)];
  Gr[l_, r_, s_] := If[phys, -(s^2/(2 l + 1)) Piecewise[{{s^l/r^(l + 1), s < r}}, r^l/s^(l + 1)],       (* free space *)
     s^2/(2 l + 1) (r^l s^(-l - 1) - s^l r^(-l - 1))];                                         (* zero Cauchy data at r = 1 *)
  cRowsC = {};
  Do[Module[{fRe = Function @@ {{zVar}, N[sRealE[l], 50]}, scale = scaleA[l], Fx, dFx, t},
     (* exact G0 = Im F~(r0+i0), G1 = Im F~'(r0+i0) of the exact F, by principal values (for the error budget only) *)
     Fx = Function @@ {{zVar}, N[(zVar - r0)^2 ssSourceModeExpr["Dzz", l, zVar, "real"], 60]};
     dFx = Function @@ {{zVar}, N[D[(zVar - r0)^2 ssSourceModeExpr["Dzz", l, zVar, "real"], zVar], 60]};
     gEx[l] = Quiet[N[{-(1/Pi) NIntegrate[Fx[t]/(t - r0), {t, 4, r0, 6}, Method -> "PrincipalValue", WorkingPrecision -> 40, PrecisionGoal -> 20],
         -(1/Pi) (Fx[4]/(4 - r0) - Fx[6]/(6 - r0) + NIntegrate[dFx[t]/(t - r0), {t, 4, r0, 6}, Method -> "PrincipalValue", WorkingPrecision -> 40, PrecisionGoal -> 20])}], NIntegrate::ncvb];      (* ncvb: ln|t - r0| in F' slows the PV; value good to < 1e-10 *)
     Do[Module[{tcs = 4 + (Range[2/h] - 1/2) h, F0, xi},
        F0 = N[(# - r0)^2 fRe[#] & /@ tcs]; SeedRandom[1000 l + Round[1/h]]; xi = RandomVariate[NormalDistribution[], Length[tcs]];
        Do[Module[{F = F0 (1 + eps xi), cs, csD, Fend, G0, G1, sF, sP, sU, eP, eU, pred, eX},
           cs = pieceSetup[F, h, 7]; csD = derivSetup[cs];
           Fend = {cs["coef"][[1]] . (-1)^Range[0, 7], Total[cs["coef"][[-1]]]};                 (* F(4), F(6) *)
           G0 = Im[cauchyAtP[cs, r0 + 0. I]];
           G1 = Im[(Fend[[1]]/(4 - r0 - I 10^-30) - Fend[[2]]/(6 - r0 - I 10^-30))/(Pi I) + cauchyAtP[csD, r0 + 0. I]];
           sF = cauchyAtP[cs, #] & /@ zsM;
           sP = (sF - I G0 - I G1 (zsM - r0))/(zsM - r0)^2; sU = sF/(zsM - r0)^2;
           eP = extErr[l, sP - sCfM[l]]; eU = extErr[l, sU - sCfM[l]];
           eX = extErr[l, (sF - I gEx[l][[1]] - I gEx[l][[2]] (zsM - r0))/(zsM - r0)^2 - sCfM[l]];       (* pinned with the exact G0, G1 *)
           pred = (Pi G1 Gr[l, #, r0] + Pi G0 (D[Gr[l, #, s], s] /. s -> r0)) & /@ If[phys, gCom, gOut];      (* -Pi G0 delta' + Pi G1 delta, as responses *)
           AppendTo[cRowsC, {l, h, eps, Max[Abs[eP]]/scale, Max[Abs[eU]]/scale, Abs[G0 - gEx[l][[1]]], Abs[G1 - gEx[l][[2]]], Max[Abs[eU - eP - pred]]/Max[Abs[pred]], Max[Abs[eX]]/scale}]],
         {eps, {0, 10^-12, 10^-8}}]], {h, hList}];
     Print["    l=", l, " [", tick[], " s]"]], {l, lA}];
  mxC[h_, eps_, k_] := Max[Cases[cRowsC, {_, h, eps, a__} :> {a}[[k]]]];
  csCheck["(C) pinned with the exact G0, G1: error falls at least as h^2.5 (shortfall in the 1/100 -> 1/200 order)", Max[0, 5/2 - Log[2, mxC[1/100, 0, 6]/mxC[1/200, 0, 6]]], 0];
  csCheck["(C) pinned with the data's G0, G1: error is first order in h (shortfall of the 1/100 -> 1/200 order below 0.8)", Max[0, 0.8 - Log[2, mxC[1/100, 0, 1]/mxC[1/200, 0, 1]]], 0];
  csCheck["(C) unpinned minus pinned = response to -Pi G0 delta' + Pi G1 delta (max relative deviation)", Max[cRowsC[[All, 8]]], 10^-6];
  mdC = StringRiffle[Join[{"| h | eps | pinned (data G0, G1) | unpinned | pinned (exact G0, G1) | error of G0 | error of G1 | noise part (pinned) |", "|---|---|---|---|---|---|---|---|"},
      Flatten@Table["| " <> ToString[h, InputForm] <> " | " <> sci[N[eps]] <> " | " <> sci[mxC[h, eps, 1]] <> " | " <> sci[mxC[h, eps, 2]] <> " | " <> sci[mxC[h, eps, 6]] <> " | " <> sci[mxC[h, eps, 3]] <> " | " <> sci[mxC[h, eps, 4]] <> " | " <>
         If[eps == 0, "-", sci[Max[Table[Abs[Cases[cRowsC, {l, h, eps, a_, __} :> a][[1]] - Cases[cRowsC, {l, h, 0, a_, __} :> a][[1]]], {l, lA}]]]] <> " |", {h, hList}, {eps, {0, 10^-12, 10^-8}}]], "\n"];
  Print["\n", mdC];
  Export[FileNameJoin[{resDir, tagR <> "C_radial_factor.csv"}], Prepend[{#[[1]], ToString[#[[2]], InputForm], sci[N[#[[3]]]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]], sci[#[[9]]]} & /@ cRowsC,
     {"l", "h", "eps", "err_pinned", "err_unpinned", "G0_error", "G1_error", "unpinned_prediction_rel_dev", "err_pinned_exactG"}], "CSV"];
  Print["    (C) done [", tick[], " s]"]]];

(* combined SS figure: (a) noise-free error against h; (b) noise-induced error against h at eps = 1e-8 *)
If[doA && doC, Module[{st = <|"A" -> Directive[Gray, Dashed], "AB" -> Directive[RGBColor[0.85, 0.45, 0.1]], "ABC" -> Directive[RGBColor[0.1, 0.35, 0.75]], "ABC+" -> Directive[RGBColor[0.15, 0.6, 0.25]]|>,
    cC = RGBColor[0.75, 0.1, 0.45], p1, p2, noiseA, noiseC},
   noiseA[v_, model_, h_, eps_] := Max[Cases[nRows, {_, v, h, model, eps, _, b_} :> b]];
   noiseC[h_, eps_] := Max[Table[Abs[Cases[cRowsC, {l, h, eps, a_, __} :> a][[1]] - Cases[cRowsC, {l, h, 0, a_, __} :> a][[1]]], {l, lA}]];
   p1 = ListLogLogPlot[Join[Table[Table[{h, mxA[v, h, 7, 1]}, {h, hList}], {v, vNames}], {Table[{h, mxC[h, 0, 1]}, {h, hList}]}], Joined -> True, PlotMarkers -> {"\[FilledCircle]", 5},
     PlotStyle -> Join[Directive[#, AbsoluteThickness[1.2]] & /@ Values[st], {Directive[cC, AbsoluteThickness[1.2]]}],
     Frame -> True, Axes -> False, FrameLabel -> {"grid spacing h", "scaled exterior SS error (no noise)"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     PlotLegends -> Placed[LineLegend[Join[Values[st], {cC}], {"subtract {A}", "subtract {A,B}", "subtract {A,B,C\:2080}", "subtract {A,B,C\:2080,C\:2081,C\:2082}", "(C) (r\[Minus]r\:2080)\.b2SS, pinned"}, LabelStyle -> $fs], Below]];
   p2 = ListLogLogPlot[{Table[{h, noiseA["ABC+", "rel. full SS", h, 10^-8]}, {h, hList}], Table[{h, noiseA["ABC+", "rel. remainder", h, 10^-8]}, {h, hList}], Table[{h, noiseC[h, 10^-8]}, {h, hList}],
      Table[{h, noiseA["ABC+", "rel. full SS", h, 10^-12]}, {h, hList}], Table[{h, noiseC[h, 10^-12]}, {h, hList}]}, Joined -> True, PlotMarkers -> {"\[FilledCircle]", 5},
     PlotStyle -> {Directive[RGBColor[0.15, 0.6, 0.25], AbsoluteThickness[1.2]], Directive[RGBColor[0.15, 0.6, 0.25], Dashed, AbsoluteThickness[1.2]], Directive[cC, AbsoluteThickness[1.2]],
       Directive[RGBColor[0.15, 0.6, 0.25], Dotted, AbsoluteThickness[1.2]], Directive[cC, Dotted, AbsoluteThickness[1.2]]},
     Frame -> True, Axes -> False, FrameLabel -> {"grid spacing h", "noise-induced exterior error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     PlotLegends -> Placed[LineLegend[{Directive[RGBColor[0.15, 0.6, 0.25]], Directive[RGBColor[0.15, 0.6, 0.25], Dashed], cC, Directive[RGBColor[0.15, 0.6, 0.25], Dotted], Directive[cC, Dotted]},
        {"subtraction, \[Epsilon] = 10\:207b\:2078 \[Times] SS", "subtraction, \[Epsilon] = 10\:207b\:2078 \[Times] remainder", "(C), \[Epsilon] = 10\:207b\:2078 \[Times] SS", "subtraction, \[Epsilon] = 10\:207b\.b9\.b2 \[Times] SS", "(C), \[Epsilon] = 10\:207b\.b9\.b2 \[Times] SS"}, LabelStyle -> $fs], Below]];
   Export[FileNameJoin[{figDir, If[phys, "fig_ss_reuse.pdf", "fig_ss_reuse_cauchy.pdf"]}], Grid[{{p1, p2}}, Spacings -> 1, Alignment -> Top]];
   Export[FileNameJoin[{figDir, If[phys, "fig_ss_reuse.csv", "fig_ss_reuse_cauchy.csv"]}], Join[{{"variant", "h", "eps", "noise_model", "max_l_scaled_error"}},
      Flatten[Table[{"subtract " <> v, ToString[h, InputForm], 0, "-", N[mxA[v, h, 7, 1], 6]}, {v, vNames}, {h, hList}], 1],
      Flatten[Table[{"radial factor (C)", ToString[h, InputForm], N[eps], "rel. full SS", N[mxC[h, eps, 1], 6]}, {h, hList}, {eps, {0, 10^-12, 10^-8}}], 1],
      Flatten[Table[{"subtract " <> v, ToString[h, InputForm], N[eps], model, N[mxN[v, h, model, eps, 1], 6]}, {v, {"ABC", "ABC+"}}, {h, hList}, {model, {"rel. remainder", "rel. full SS"}}, {eps, epsList}], 3]], "CSV"]]];

(* ======================================================================================= *)
(* (D) complex-r angular quadrature (no subtraction, no Cauchy step)                        *)
(* ======================================================================================= *)
If[doD, Module[{},
  Print["\n--- (D) complex-r angular quadrature ---"];
  (* mu -> rho = |x - x0| (rho^2 = z^2 + r0^2 - 2 z r0 mu, dmu = -rho drho/(z r0)): the integrand is rational in rho with its only pole at rho = 0.
     The mu path from -1 to 1 becomes a rho path from z + r0 to r0 - z (in) or z - r0 (out); via -i/2 (in) or +i/2 (out) it never meets rho = 0,
     so the result is the continuation from a_L (in) or from a_R (out).  The deformation of the theta path is implicit in this map. *)
  nq = 80; {tq, wq} = Transpose[N[GaussianQuadratureWeights[nq, 0, 1, 30]]];
  seg[g_, a_, b_] := (b - a) Total[wq (g /@ (a + (b - a) tq))];
  pathInt[g_, z_, br_] := With[{via = If[br === "in", -I/2, I/2], e = If[br === "in", r0 - z, z - r0]}, seg[g, z + r0, via] + seg[g, via, e]];
  muOf[z_, rho_] := (z^2 + r0^2 - rho^2)/(2 z r0); dOf[z_, rho_] := (z^2 - r0^2 - rho^2)/(2 r0);          (* d = z mu - r0 *)
  gSS[l_, z_] := Function[rho, (1/2) q^2 (-2/rho^4 + 8 dOf[z, rho]^2/rho^6) LegendreP[l, muOf[z, rho]] (-rho/(z r0))];
  gSR[l_, z_] := Function[rho, q ((c0 + c1 z muOf[z, rho]) (3 dOf[z, rho]^2 - rho^2)/rho^5 - 2 c1 dOf[z, rho]/rho^3) LegendreP[l, muOf[z, rho]] (-rho/(z r0))];
  modeD[g_, l_, z_, br_] := Sqrt[Pi (2 l + 1)] pathInt[g[l, z], z, br];
  Hm[z_] := 1 + (I/Pi) Log[z - r0];
  (* undifferentiated product modes at complex z with r-derivatives at fixed mu, along the rho path of branch br *)
  pathPts[z_, br_] := With[{via = If[br === "in", -I/2, I/2], e = If[br === "in", r0 - z, z - r0]},
     Join[Transpose[{z + r0 + (via - z - r0) tq, (via - z - r0) wq}], Transpose[{via + (e - via) tq, (e - via) wq}]]];
  prodModes[z_, br_, h3_, jm_] := Module[{pts = pathPts[z, br], rho, w, mu, vals, leg},
     rho = pts[[All, 1]]; w = pts[[All, 2]]; mu = muOf[z, rho];
     vals = (# (-rho/(z r0)) w) & /@ h3[z, mu, rho];
     leg = Table[LegendreP[j, mu], {j, 0, jm}];
     Table[Sqrt[Pi (2 j + 1)] (leg[[j + 1]] . #) & /@ vals, {j, 0, jm}]];
  (* SR = (1/2) D_z^2 of [2 q PhiR/rho]_j; the Theta -> H blend in + (out - in) H acts on the product, before D_z^2 *)
  hSR[z_, mu_, rho_] := Module[{u = z - r0 mu, PR = c0 + c1 z mu}, 2 q {PR/rho, c1 mu/rho - PR u/rho^3, -2 c1 mu u/rho^3 + PR (-1/rho^3 + 3 u^2/rho^5)}];
  blendDzz[l_, z_, h3_] := Module[{pin = prodModes[z, "in", h3, l + 2], pout = prodModes[z, "out", h3, l + 2], H0 = Hm[z], H1 = I/(Pi (z - r0)), H2 = -I/(Pi (z - r0)^2), F},
     F = Table[Module[{a = pin[[j + 1]], d = pout[[j + 1]] - pin[[j + 1]]}, {a[[1]] + d[[1]] H0, a[[2]] + d[[2]] H0 + d[[1]] H1, a[[3]] + d[[3]] H0 + 2 d[[2]] H1 + d[[1]] H2}], {j, 0, l + 2}];
     (1/2) dzzArr[F, z][[l + 1, 1]]];
  dRows = {};
  Do[Module[{fSRc = Function @@ {{zVar}, N[(1/2) dzzExpr[Function[j, If[j < 0, 0, 2 q Sum[KCouplingContinued[j, lp, zVar, "ThetaToH"] phiRAnalytic[lp][zVar], {lp, 0, 1}]]], l], 40]},
       ssIn, ssOut, srIn, srOut, srBl, srCf, scSR, eSSin, eSSout, eSRin, eSRbl, srPB, sub = Range[1, Length[zsM], 4]},
     ssIn = modeD[gSS, l, #, "in"] & /@ zsM; ssOut = modeD[gSS, l, #, "out"] & /@ zsM;
     srIn = modeD[gSR, l, #, "in"] & /@ zsM; srOut = modeD[gSR, l, #, "out"] & /@ zsM; srBl = srIn + (srOut - srIn) (Hm /@ zsM);
     srCf = N[fSRc /@ zsM];
     srPB = blendDzz[l, #, hSR] & /@ zsM[[sub]];
     scSR = Max[Abs[analyticEval[pieceAnalyticSolution["Dzz", "SR", l], #] & /@ gCom]];
     eSSin = Max[Abs[extErr[l, ssIn - sCfM[l]]]]/scaleA[l]; eSSout = Max[Abs[extErr[l, ssOut - sCfM[l]]]]/scaleA[l];
     eSRin = Max[Abs[extErr[l, srIn - srCf]]]/scSR; eSRbl = Max[Abs[extErr[l, srBl - srCf]]]/scSR;
     AppendTo[dRows, {l, Max[Abs[ssIn - sCfM[l]]]/Max[Abs[sCfM[l]]], Max[Abs[ssOut - sCfM[l]]]/Max[Abs[sCfM[l]]], eSSin, eSSout, eSRin, eSRbl, Max[Abs[srPB - srCf[[sub]]]]/Max[Abs[srCf[[sub]]]]}];
     Print["    l=", l, " [", tick[], " s]"]], {l, lA}];
  csCheck["(D) SS by complex-r angular quadrature (continued from a_L) = closed form, exterior (max over l, scaled)", Max[dRows[[All, 4]]], 10^-11];
  csCheck["(D) SS from a_R (out branch) = closed form pointwise on the contour (max over l, relative)", Max[dRows[[All, 3]]], 10^-11];
  csCheck["(D) SR: Theta->H blend of the angular-quadrature product modes, then D_z^2 = product route, pointwise (max over l, relative)", Max[dRows[[All, 8]]], 10^-10];
  csCheck["(D) negative control: blending the differentiated SR (in, out) with H fails (exterior, max over l, scaled > 1e-2)", If[Max[dRows[[All, 7]]] > 10^-2, 0, 1], 0];
  (* cost: adaptive mu-quadrature evaluations for 1e-12, real r near r0 against contour points *)
  costRows = Module[{l = 2, cnt, val, f},
     f[z_, mu_] := (1/2) q^2 (-2/(z^2 + r0^2 - 2 z r0 mu)^2 + 8 (z mu - r0)^2/(z^2 + r0^2 - 2 z r0 mu)^3) LegendreP[l, mu];
     Join[Table[cnt = 0; val = Sqrt[Pi (2 l + 1)] NIntegrate[f[r0 - dl, mu], {mu, -1, 1}, PrecisionGoal -> 12, MaxRecursion -> 60, EvaluationMonitor :> cnt++];
        {"real r = r0 - " <> ToString[dl, InputForm], cnt, Abs[val - N[sRealE[l] /. zVar -> r0 - dl, 30]]/Abs[val]}, {dl, {1/10, 1/100, 1/400, 1/1000}}],
      Table[cnt = 0; val = Sqrt[Pi (2 l + 1)] NIntegrate[gSS[l, zz][rho], {rho, zz + r0, -I/2, r0 - zz}, PrecisionGoal -> 12, MaxRecursion -> 60, EvaluationMonitor :> cnt++];
        {"contour z = " <> ToString[N[zz, 3]], cnt, Abs[val - modeD[gSS, l, zz, "in"]]/Abs[val]}, {zz, N[yA /@ {4.6, 5, 5.4}]}]]];
  (* SR Part 1 of the hybrid split for the Stage 6 puncture (n = 3, l'max = 20, k = 2): closed form against in + (out - in) H from angular quadrature *)
  p1Row = Module[{n = 3, cc = r0^-4, M = 20, k = 2, tc, tp, A0, A1, A2, B2, gP1, l = 2, cfP1, dP1},
     tc = taylor9[n, cc, M]; tp = tPoly9[k, tc]; {A0, A1, A2, B2} = N[tc[[1 ;; 4]], 30];
     hP1[z_, mu_, rho_] := Module[{u = z - r0 mu, d = z mu - r0, s2 = z^2 (1 - mu^2), T, dT, ddT, P, dP, ddP},
        T = A0 + A1 d + A2/2 d^2 + B2/2 s2; dT = A1 mu + A2 d mu + B2 z (1 - mu^2); ddT = A2 mu^2 + B2 (1 - mu^2);
        P = q/rho + cc q rho^3; dP = -q u/rho^3 + 3 cc q rho u; ddP = q (-1/rho^3 + 3 u^2/rho^5) + 3 cc q (u^2/rho + rho);
        2 {P T, dP T + P dT, ddP T + 2 dP dT + P ddT}];
     cfP1 = Function @@ {{zVar}, N[(1/2) dzzExpr[Function[j, If[j < 0, 0, 2 Sum[KPc[n, cc, j, lp] tp[[lp + 1, 1]], {lp, 0, k}]]], l], 40]};
     dP1 = Table[Module[{zz = zsM[[i]], v}, v = N[cfP1[zz]]; Abs[blendDzz[l, zz, hP1] - v]/Abs[v]], {i, 1, Length[zsM], 7}];
     {Max[dP1]}];
  csCheck["(D) SR Part 1 (Stage 6 puncture, n = 3, k = 2, l = 2): Theta->H blend of angular-quadrature product modes, then D_z^2 = closed form (relative)", p1Row[[1]], 10^-9];
  (* continued residual modes of the Stage 6 puncture at complex z: geometric growth, so the 3D mode sum diverges *)
  divRows = Module[{n = 3, cc = r0^-4, Lm = 80}, Table[Module[{zz = N[yA[aa], 30], v, slope, ps0, psR},
       v = Table[evR[resC[n, cc, lp], zz], {lp, 0, Lm}];
       slope = Exp[Coefficient[Fit[Table[{lp, Log[Abs[v[[lp + 1]]]]}, {lp, 60, Lm}], {1, x}, x], x]];
       ps0 = Abs[Accumulate[v Table[Sqrt[(2 lp + 1)/(4 Pi)], {lp, 0, Lm}]]];
       psR = Abs[Accumulate[Table[pwVal[resR[n, cc, lp], N[Re[zz], 30]], {lp, 0, Lm}] Table[Sqrt[(2 lp + 1)/(4 Pi)], {lp, 0, Lm}]]];
       {aa, N[Abs[zz]/r0, 4], N[slope, 4], ps0[[{21, 41, 61, 81}]], psR[[{21, 41, 61, 81}]], Abs[v[[Lm + 1]]], Ordering[Abs[v[[2 ;;]]], 1][[1]]}], {aa, {4.6, 5., 5.4, 5.49}}]];
  mdD = StringRiffle[Join[{"| l | SS in: pointwise rel. diff | SS out: pointwise | SS in: exterior | SS out: exterior | SR in only: exterior | SR, H blend after D_z^2: exterior | SR, H blend of products then D_z^2: pointwise |", "|---|---|---|---|---|---|---|---|"},
      ("| " <> ToString[#[[1]]] <> " | " <> StringRiffle[sci /@ #[[2 ;;]], " | "] <> " |") & /@ dRows], "\n"];
  mdD2 = StringRiffle[Join[{"| point | adaptive evaluations for 1e-12 | relative error |", "|---|---|---|"}, ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> sci[#[[3]]] <> " |") & /@ costRows], "\n"];
  mdD3 = StringRiffle[Join[{"| a | abs(z)/r0 | growth per l' (fit 60-80) | smallest term at l' | abs(term) at l' = 80 | abs(partial sum at theta=0), l' = 20, 40, 60, 80 | real r = Re z: same |", "|---|---|---|---|---|---|---|"},
      ("| " <> ToString[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> ToString[#[[7]]] <> " | " <> sci[#[[6]]] <> " | " <> StringRiffle[sci /@ #[[4]], ", "] <> " | " <> StringRiffle[sci /@ #[[5]], ", "] <> " |") & /@ divRows], "\n"];
  Print["\n", mdD, "\n\n", mdD2, "\n\nSR Part 1 (n = 3, k = 2, l = 2), blend vs closed form: ", sci[p1Row[[1]]], "\n\n", mdD3];
  Export[FileNameJoin[{resDir, tagR <> "D_angular.csv"}], Prepend[{#[[1]], Sequence @@ (sci /@ #[[2 ;;]])} & /@ dRows, {"l", "SS_in_pointwise", "SS_out_pointwise", "SS_in_exterior", "SS_out_exterior", "SR_in_exterior", "SR_blend_after_Dzz_exterior", "SR_product_blend_pointwise"}], "CSV"];
  Export[FileNameJoin[{resDir, tagR <> "D_summary.md"}], "# Stage 12 (D): complex-r angular quadrature\n\n" <> mdD <> "\n\n## cost (l = 2)\n\n" <> mdD2 <>
     "\n\nSR Part 1 (Stage 6 puncture, n = 3, k = 2, l = 2): in + (out - in) H from angular quadrature vs closed form, max relative difference " <> sci[p1Row[[1]]] <>
     "\n\n## continued residual modes of the Stage 6 puncture (n = 3), summed at complex z\n\n" <> mdD3 <> "\n", "Text"];
  Print["    (D) done [", tick[], " s]"]]];

(* ======================================================================================= *)
(* (B) H1 (k = 0) convergence, Stage 6 puncture                                            *)
(* ======================================================================================= *)
If[doB, Module[{},
  Print["\n--- (B) H1 (k = 0) convergence ---"];
  nB = {1, 3, 5, 7}; jMax = 8; Mm = 40; lB = Range[0, 6];
  gridPlot = Join[Subdivide[1, 9/2, 14], Subdivide[11/2, 10, 18]]; gAll = Union[gCom, gridPlot];
  posCom = Flatten[Position[gAll, #] & /@ gCom]; posOut = Flatten[Position[gAll, #] & /@ gOut];
  nodesB = cNodes[12, gAll];
  (* cumulative (over l') F_j for SR and RR at one point *)
  nodeFcum[{A_, V_}, Mx_] := Module[{srC, rrSh, tr, z3 = 0 V[[1]]},
     srC = Table[2 prodD[Sum[gt[j, lp, L] A[[L + 1]], {L, Abs[j - lp], j + lp}], V[[lp + 1]]], {j, 0, jMax}, {lp, 0, Mx}];
     tr[j_, l1_, l2_] := With[{g = gt[j, l1, l2]}, If[g == 0, z3, g prodD[V[[l1 + 1]], V[[l2 + 1]]]]];
     rrSh = Table[Sum[tr[j, m, l2], {l2, Max[0, m - j], m}] + Sum[tr[j, l1, m], {l1, Max[0, m - j], m - 1}], {j, 0, jMax}, {m, 0, Mx}];
     {Accumulate /@ srC, Accumulate /@ rrSh}];
  LamU[l_, nn_] := 2^nn (2 l + 1)/((2 l + nn + 1) Pochhammer[l + nn/2 + 3/2, -nn]);          (* Upton et al. eq. (173), s = 0 *)
  csCheck["(B) Upton Lambda_{l,0,-4}: sum over l >= 0 vanishes (eq. 174)", Abs[N[Sum[LamU[l, -4], {l, 0, Infinity}], 30]], 10^-25];
  eB = <||>; dB = <||>; slopes = {}; accRows = {};
  Do[Module[{c = r0^-(n + 1), P, AVin, AVout, cin, cout, a0, ep},
     P = Table[Module[{cum = nodeFcum[nodeAV[n, c, s, Automatic, jMax + Mm, Mm], Mm]},
         Table[(1/2) dzzArr[cum[[1, All, M + 1]] + cum[[2, All, M + 1]], s][[All, 1]], {M, 0, Mm}]], {s, nodesB[[All, 1]]}];
     AVin = nodeAV[n, c, r0, "in", jMax + Mm, Mm, 4]; AVout = nodeAV[n, c, r0, "out", jMax + Mm, Mm, 4];
     cin = nodeFcum[AVin, Mm]; cout = nodeFcum[AVout, Mm];
     a0 = Accumulate[Table[Y0[l] evR[resR[n, c, l][[1]], r0], {l, 0, Mm}]];                  (* T_0 = Phi^res_M(x0), mode sum *)
     ep = exPairF[n, c];
     Print["    n=", n, ": ", Length[nodesB], " nodes cached [", tick[], " s]"];
     Do[Module[{exE = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, ep[j][[b]]]], l], {b, 2}], exPt, gEx, g0, scale},
        exPt = pwVal[exE, #] & /@ nodesB[[All, 1]]; gEx = (1/2) N[dzzDeltaCoeff[ep, l], 40]; g0 = G0[l, #, r0] & /@ gAll;
        scale = Max[Abs[analyticEval[Sum[pieceAnalyticSolution["Dzz", p, l], {p, {"SS", "SR", "RR"}}], #] & /@ gCom]];
        (* H1, k = 0: real-line distributional solve minus the delta Part 2 misses: e = response(trunc - exact) + (gP1 + gRR - gEx) G0 *)
        eB[{n, l}] = Table[Module[{gRR = (1/2) dzzDeltaNum[l, cin[[2, # + 1, M + 1, 2]] &, cout[[2, # + 1, M + 1, 2]] &], tp = tPoly9[0, {a0[[M + 1]], 0, 0, 0, 0, 0}], gP1, gSR},
            gP1 = (1/2) dzzDeltaNum[l, p1F[AVin[[1]], tp, 0, r0][[# + 1, 2]] &, p1F[AVout[[1]], tp, 0, r0][[# + 1, 2]] &];
            gSR = (1/2) dzzDeltaNum[l, cin[[1, # + 1, M + 1, 2]] &, cout[[1, # + 1, M + 1, 2]] &];
            If[M == Mm, csCheck["(B) n=" <> ToString[n] <> " l=" <> ToString[l] <> ": k = 0 Part 1 carries the whole SR delta (relative)", Abs[gSR - gP1]/Abs[gSR], 10^-25]];
            cResponse[l, nodesB, P[[All, M + 1, l + 1]] - exPt, gAll] + (gP1 + gRR - gEx) g0], {M, 0, Mm}]/scale;
        (* per-l' contributions: d_l' = e(l') - e(l'-1); e(-1) = minus the exact field *)
        dB[{n, l}] = Prepend[Differences[eB[{n, l}]], eB[{n, l}][[1]] - (cResponse[l, nodesB, -exPt, gAll] - gEx g0)/scale]], {l, lB}];
     Print["    n=", n, " done [", tick[], " s]"]], {n, nB}];
  errM[n_, M_] := Max[Table[Max[Abs[eB[{n, l}][[M + 1, posCom]]]], {l, lB}]];
  (* slopes: least squares of log(max error) on log(l'max) over l'max in [20, 40]; per-mode power from |d_l'| at r = 6.5 *)
  Do[Module[{fitE = Fit[Table[{Log[M], Log[errM[n, M]]}, {M, 20, Mm}], {1, x}, x], lmx, fitD},
     lmx = First@Ordering[Table[-Abs[dB[{n, l}][[Mm + 1, posOut[[1]]]]], {l, lB}], 1] - 1;
     fitD = Fit[Table[{Log[lp], Log[Abs[dB[{n, lmx}][[lp + 1, posOut[[1]]]]]]}, {lp, 25, Mm}], {1, x}, x];
     AppendTo[slopes, {n, -Coefficient[fitE, x], -Coefficient[fitD, x], lmx, errM[n, 12], errM[n, 20], errM[n, 40]}]], {n, nB}];
  Do[csCheck["(B) n=" <> ToString[s[[1]]] <> ": error slope = n + 2 within 0.3 (fit over l'max 20..40)", Abs[s[[2]] - (s[[1]] + 2)], 3/10], {s, slopes}];
  Do[csCheck["(B) n=" <> ToString[s[[1]]] <> ": per-mode power = n + 3 within 0.3 (r = 6.5, l'max 25..40)", Abs[s[[3]] - (s[[1]] + 3)], 3/10], {s, slopes}];
  (* tail acceleration at each l'max: fit d_l' on l' in [M-5, M] with 2 terms of known leading power p = n + 3, add the tail analytically *)
  accel[n_, l_, M_, basis_] := Module[{pw = n + 3, bf, X, Y, coef, tl, wn = If[basis === "Lambda1", 3, 6]},
     bf = Which[basis === "power", Function[lp, {lp^-pw, lp^(-pw - 1)}], basis === "Lambda", Function[lp, {LamU[lp, -pw], LamU[lp, -pw - 2]}], basis === "Lambda1", Function[lp, {LamU[lp, -pw]}]];
     X = N[bf /@ Range[M - wn + 1, M], 30]; Y = dB[{n, l}][[M - wn + 2 ;; M + 1]];
     coef = PseudoInverse[X] . Y;                                                           (* per exterior point *)
     tl = If[basis === "power", N[{HurwitzZeta[pw, M + 1], HurwitzZeta[pw + 1, M + 1]}, 30], -N[Sum[bf[lp], {lp, 0, M}], 30]];   (* sum_{l'>M} Lambda = -sum_{l'<=M} *)
     eB[{n, l}][[M + 1]] + tl . coef];
  accM[n_, M_, basis_] := Max[Table[Max[Abs[accel[n, l, M, basis][[posCom]]]], {l, lB}]];
  Do[AppendTo[accRows, {n, M, errM[n, M], accM[n, M, "power"], accM[n, M, "Lambda"], accM[n, M, "Lambda1"]}], {n, nB}, {M, 8, Mm}];
  (* validation by a contour solve for n = 1 (new here): H1 k = 0, l'max = 8, l = 2 *)
  Module[{n = 1, c = r0^-2, M = 8, l = 2, LA, tc, tCh, aN, yAc, srAt, sCh, sG, sL, AVin, AVout, FIn, FOut, gSR, lineAR, tp, p1At, p1Ch, p1G, src, gP1, d, pred},
   LA = jMax + M; tc = taylor9[n, c, M];
   srAt[s_] := pt[l, nodeF[nodeAV[n, c, s, Automatic, LA, M], {M}, False][M][[1]], s];
   AVin = nodeAV[n, c, r0, "in", LA, M]; AVout = nodeAV[n, c, r0, "out", LA, M];
   FIn = nodeF[AVin, {M}, False][M][[1]]; FOut = nodeF[AVout, {M}, False][M][[1]];
   gSR = (1/2) dzzDeltaNum[l, FIn[[# + 1, 2]] &, FOut[[# + 1, 2]] &]; lineAR = gSR {G0[l, aR, r0], dG0[l, aR, r0]};
   tCh = {chebX[cI[[1]], r0, NC], chebX[r0, cI[[2]], NC]}; aN = aCheb[NA]; yAc = N[yC[#], 50] & /@ aN;
   sCh = Map[srAt, tCh, {2}]; sG = Map[srAt, gaussP[[All, All, 1]], {2}]; sL = srAt /@ lineNodes[[All, 1]];
   tp = tPoly9[0, tc]; p1At[s_] := pt[l, p1F[nodeAV[n, c, s, Automatic, LA, 0][[1]], tp, 0, s], s];
   p1Ch = Map[p1At, tCh, {2}]; p1G = Map[p1At, gaussP[[All, All, 1]], {2}];
   src = (evR[(1/2) dzzExpr[Function[j, If[j < 0, 0, 2 KPc[n, c, j, 0] tp[[1, 1]]]], l], #] & /@ yAc) + (cauchyAt[sG - p1G, derivCoefs /@ (sCh - p1Ch), #] & /@ yAc);
   gP1 = (1/2) dzzDeltaNum[l, p1F[AVin[[1]], tp, 0, r0][[# + 1, 2]] &, p1F[AVout[[1]], tp, 0, r0][[# + 1, 2]] &];
   d = deltaAt[l, src, sL] - lineAR; pred = -(gSR - gP1) {G0[l, aR, r0], dG0[l, aR, r0]};
   valB = {d, pred, Max[Abs[lineAR]]};
   csCheck["(B) validation n=1, l'max=8, l=2: H1 k=0 contour - distributional line at aR = 0 (relative to |gSR G0|)", Max[Abs[d - pred]]/Max[Abs[lineAR]], 10^-8]];
  (* tables *)
  mdB1 = StringRiffle[Join[{"| n | slope of error (l'max 20-40) | expected n+2 | per-mode power (r = 6.5, l' 25-40) | expected n+3 | error at l'max = 12 | 20 | 40 |", "|---|---|---|---|---|---|---|---|"},
      ("| " <> ToString[#[[1]]] <> " | " <> ToString[NumberForm[#[[2]], 3]] <> " | " <> ToString[#[[1]] + 2] <> " | " <> ToString[NumberForm[#[[3]], 3]] <> " | " <> ToString[#[[1]] + 3] <> " | " <>
         StringRiffle[sci /@ #[[5 ;; 7]], " | "] <> " |") & /@ slopes], "\n"];
  mdB2 = StringRiffle[Join[{"| n | l'max | plain | + power-law tail (2 terms) | + Lambda tail (2 terms) | + Lambda tail (1 term) | gain (power) | gain (Lambda, 2) | gain (Lambda, 1) |", "|---|---|---|---|---|---|---|---|---|"},
      ("| " <> ToString[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;; 6]], " | "] <> " | " <> StringRiffle[ToString[NumberForm[#, 3]] & /@ (#[[3]]/#[[4 ;; 6]]), " | "] <> " |") & /@
       Select[accRows, MemberQ[{12, 20, 30}, #[[2]]] &]], "\n"];
  mdB3 = StringRiffle[Join[{"| n | l'max | max error, r <= 3.5 | r in [3.5, 4.5] | r in [5.5, 6.5] | r in [6.5, 10] |", "|---|---|---|---|---|---|"},
      Flatten@Table[Module[{prof = Table[Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]], {k, Length[gAll]}], regMax},
         regMax[{lo_, hi_}] := Max[Table[If[lo <= gAll[[k]] <= hi, prof[[k]], 0], {k, Length[gAll]}]];
         "| " <> ToString[n] <> " | " <> ToString[M] <> " | " <> StringRiffle[sci[regMax[#]] & /@ {{1, 7/2}, {7/2, 9/2}, {11/2, 13/2}, {13/2, 10}}, " | "] <> " |"],
        {n, {3, 5}}, {M, {12, 20}}]], "\n"];
  Print["\n", mdB1, "\n\n", mdB2, "\n\n", mdB3, "\n\nvalidation (n = 1): contour - line = ", sci /@ valB[[1]], ", predicted ", sci /@ valB[[2]]];
  Export[FileNameJoin[{resDir, "stage12_B_summary.md"}], "# Stage 12 (B): H1 (k = 0) convergence, Stage 6 puncture, Dzz\n\nScaled exterior error, max over l = 0..6 and r in the exterior grid.\n\n" <> mdB1 <>
     "\n\n## tail acceleration (fit of d_l' on l' in [l'max-5, l'max], leading power n+3, two terms)\n\n" <> mdB2 <> "\n\n## error against r (max over l, scaled)\n\n" <> mdB3 <>
     "\n\nvalidation by a contour solve (n = 1, l'max = 8, l = 2): contour - distributional line at aR = " <> StringRiffle[sci /@ valB[[1]], ", "] <> " (predicted 0; |gSR G0| = " <> sci[valB[[3]]] <> ").\n", "Text"];
  Export[FileNameJoin[{resDir, "stage12_convergence.csv"}], Prepend[Flatten[Table[{n, M, sci[errM[n, M]]}, {n, nB}, {M, 0, Mm}], 1], {"n", "lpmax", "max_err"}], "CSV"];
  Export[FileNameJoin[{resDir, "stage12_acceleration.csv"}], Prepend[{#[[1]], #[[2]], sci[#[[3]]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]]} & /@ accRows, {"n", "lpmax", "plain", "power_tail", "Lambda_tail", "Lambda1_tail"}], "CSV"];
  Export[FileNameJoin[{resDir, "stage12_permode.csv"}], Prepend[Flatten[Table[{n, l, lp, sci[dB[{n, l}][[lp + 1, posOut[[1]]]]]}, {n, nB}, {l, lB}, {lp, 0, Mm}], 2], {"n", "l", "lp", "d_lp_at_r6.5"}], "CSV"];
  Export[FileNameJoin[{resDir, "stage12_profile.csv"}], Prepend[Flatten[Table[{n, M, N[gAll[[k]], 6], sci[Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]]]}, {n, nB}, {M, {12, 20, 40}}, {k, Length[gAll]}], 2], {"n", "lpmax", "r", "max_l_err"}], "CSV"];
  (* figure: (a) error vs l'max, (b) per-mode contributions, (c) error vs r, (d) acceleration *)
  Module[{cols = <|1 -> RGBColor[0.55, 0.3, 0.7], 3 -> RGBColor[0.1, 0.35, 0.75], 5 -> RGBColor[0.15, 0.6, 0.25], 7 -> RGBColor[0.85, 0.45, 0.1]|>, pa, pb, pc, pd, fl},
   fl[s_, n_] := Module[{a = Log[errM[n, 30]] + s Log[30]}, Table[{M, Exp[a - s Log[M]]}, {M, {15, Mm}}]];
   pa = ListLogLogPlot[Join[Table[Table[{M, errM[n, M]}, {M, 2, Mm}], {n, nB}], Table[fl[Select[slopes, #[[1]] == n &][[1, 2]], n], {n, nB}]], Joined -> True,
     PlotStyle -> Join[Table[Directive[cols[n], AbsoluteThickness[1.3]], {n, nB}], Table[Directive[Black, Dashed, AbsoluteThickness[0.6]], {n, nB}]],
     Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]\[Prime]\:2098\:2090\:2093", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     Epilog -> Table[Inset[Style["n=" <> ToString[s[[1]]] <> ": slope " <> ToString[NumberForm[s[[2]], 3]], $fs, cols[s[[1]]]], Scaled[{0.03, 0.03 + 0.08 (4 - Position[nB, s[[1]]][[1, 1]])}], {Left, Bottom}], {s, slopes}]];
   pb = ListLogLogPlot[Table[Table[{lp, Abs[dB[{s[[1]], s[[4]]}][[lp + 1, posOut[[1]]]]]}, {lp, 2, Mm}], {s, slopes}], Joined -> True,
     PlotStyle -> Table[Directive[cols[n], AbsoluteThickness[1.3]], {n, nB}], Frame -> True, Axes -> False,
     FrameLabel -> {"\[ScriptL]\[Prime]", "|per-mode contribution| at r = 6.5"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     Epilog -> Table[Inset[Style["n=" <> ToString[s[[1]]] <> ": power " <> ToString[NumberForm[s[[3]], 3]], $fs, cols[s[[1]]]], Scaled[{0.03, 0.03 + 0.08 (4 - Position[nB, s[[1]]][[1, 1]])}], {Left, Bottom}], {s, slopes}]];
   pc = ListLogPlot[Flatten[Table[{Table[{gAll[[k]], Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]] + 10^-30}, {k, Select[Range[Length[gAll]], gAll[[#]] <= 9/2 &]}],
         Table[{gAll[[k]], Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]] + 10^-30}, {k, Select[Range[Length[gAll]], gAll[[#]] >= 11/2 &]}]}, {n, {3, 5}}, {M, {12, 20}}], 2],
     Joined -> True, PlotStyle -> Flatten[Table[{#, #} &[Directive[cols[n], If[M == 12, Dashed, Opacity[1]], AbsoluteThickness[1.2]]], {n, {3, 5}}, {M, {12, 20}}]],
     Frame -> True, Axes -> False, FrameLabel -> {"r", "max\:2097 scaled error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> {{1, 10}, All},
     Epilog -> {GrayLevel[0.85], Rectangle[Scaled[{0.39, 0}], Scaled[{0.5, 1}]], Black, Inset[Style["tube", $fs], Scaled[{0.445, 0.5}]],
       Inset[Style["dashed: \[ScriptL]\[Prime]\:2098\:2090\:2093 = 12, solid: 20", $fs], Scaled[{0.98, 0.03}], {Right, Bottom}]}];
   pd = ListLogLogPlot[Flatten[Table[{{#[[2]], #[[3]]} & /@ Select[accRows, #[[1]] == n &], {#[[2]], #[[5]]} & /@ Select[accRows, #[[1]] == n &], {#[[2]], #[[6]]} & /@ Select[accRows, #[[1]] == n &]}, {n, {3, 5}}], 1],
     Joined -> True, PlotStyle -> Flatten[Table[{Directive[cols[n], AbsoluteThickness[1.3]], Directive[cols[n], Dashed, AbsoluteThickness[1.3]], Directive[cols[n], Dotted, AbsoluteThickness[1.3]]}, {n, {3, 5}}], 1],
     Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]\[Prime]\:2098\:2090\:2093", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     Epilog -> {Inset[Style["n = 3 (blue), 5 (green): solid plain;", $fs], Scaled[{0.03, 0.11}], {Left, Bottom}], Inset[Style["dashed + \[CapitalLambda] tail, 2 terms; dotted, 1 term", $fs], Scaled[{0.03, 0.03}], {Left, Bottom}]}];
   Export[FileNameJoin[{figDir, "fig_hybrid_convergence_cauchy.pdf"}], Grid[{{pa, pb}, {pc, pd}}, Spacings -> {1, 1}]];          (* Cauchy data at r = 1; the PRL figure (physical conditions) is made by stage12b *)
   Export[FileNameJoin[{figDir, "fig_hybrid_convergence_cauchy.csv"}], Prepend[{#[[1]], #[[2]], N[#[[3]], 6], N[#[[4]], 6], N[#[[5]], 6], N[#[[6]], 6]} & /@ accRows, {"n", "lpmax", "plain", "power_tail", "Lambda_tail", "Lambda1_tail"}], "CSV"]];
  Print["    (B) done [", tick[], " s]"]]];

Print["\nStage 12 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, If[phys, "stage12phys_checks.csv", "stage12_checks" <> If[doA && doB, "", If[doA, "_ACD", "_B"]] <> ".csv"]}],
  Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
