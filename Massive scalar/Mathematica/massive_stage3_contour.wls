#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Massive scalar, Stage 3: the second-order sources on the contour and the exterior field with the exact puncture q e^{-m rho}/rho.
   (A) T2: the SS source (1/2) d_Z^2 (q^2 e^{-2 m rho}/rho^2) at the contour nodes by complex-radius angular quadrature in rho.
       The rho-integrand is e^{-2 m rho} x (Laurent polynomial): still one fixed pole at rho = 0, but with odd powers, so the two branches
       differ by a real-analytic function (kinks) as well as by I Pi x residue (log).  The modes are therefore blended, Fc + (Fout - Fc) H
       before the d_Z^2 ladder (Fc = in-branch quadrature + I Pi A, A the residue at rho = 0), which also produces the new delta term
       4 Pi m q^2 Y_l(0)/r0^2 of the massive SS source.  Checked against the closed-form continuation (exponential integrals) pointwise and in
       the exterior field (contour Green's-function quadrature with the physical conditions); node-count study; single-branch control.
   (B) the SR source 2Q[q e^{-m rho}/rho, PhiR] by the same quadrature (blend of the undifferentiated product modes) against the closed form.
   (C) exact-puncture accuracy: three-region NDSolve solves (real axis - contour - real axis, physical conditions) of SS and SR with the
       closed-form contour sources, RR on the real line; errors against the Stage 2 reference; point-supported fit (T3 (iv)) with a
       positive control; the pointwise real-line SR without its delta (negative control) and with c_l G(r, r0) restored.
   (D) noise: 1e-8 relative noise on the SS mode samples, contour nodes against a real-axis grid graded toward r0 (noise-induced error only).
   (E) m = 0 against the massless Stage 9 (Test 1) numbers.  Figures: contour and l = 2 sources (PRD Fig. 1), source modes, recovered modes.
   Run:  wolframscript -file massive_stage3_contour.wls [smoke]                                                                       *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourMassive.wl"}]];
resDir = FileNameJoin[{ParentDirectory[dir], "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
smoke = MemberQ[$ScriptCommandLine, "smoke"];
csResetChecks[]; t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
Off[N::meprec]; Off[General::munfl]; Off[NIntegrate::ncvb];
mList = If[smoke, {1}, {1/5, 1, 0}]; lmax = If[smoke, 2, 6]; lList = Range[0, lmax];
r0 = 5; q = 1; c0 = 1; c1 = 1/10; rMin = 1; rMax = 10; aL = 9/2; aR = 11/2;
gIn = Subdivide[1, 7/2, 12]; gOut = Subdivide[13/2, 10, 14]; gCom = Join[gIn, gOut];
mStr[m_] := ToString[m, InputForm]; mTag[m_] := StringReplace[mStr[m], "/" -> "_"];
$fs = {FontFamily -> "Times", FontSize -> 9};
Print["=== Massive Stage 3: contour sources and exact-puncture accuracy; m = ", mList, ", l <= ", lmax, " ==="];

(* ---------------- reference (Stage 2) ---------------- *)
loadRef[m_] := Module[{f = FileNameJoin[{resDir, "massive_reference_m" <> mTag[m] <> ".mx"}], ref, ok = False},
   If[FileExistsQ[f], ref = Import[f]; ok = Quiet[TrueQ[AssociationQ[ref] && KeyExistsQ[ref, "total"] && Max[Keys[ref["total"]]] >= lmax &&
        Abs[robinResidual[0, ref["total"][0]["val"][[1]], ref["total"][0]["der"][[1]], rMin]] < 10^-15]]];
   If[! ok, Print["    (reference for m = ", mStr[m], " missing or invalid: recomputing)"]; ref = <||>; Do[ref[p] = exactReference[lmax, gCom, {p}], {p, {"SS", "SR", "RR"}}];
    ref["total"] = Association@Table[l -> <|"val" -> Sum[ref[p][l]["val"], {p, {"SS", "SR", "RR"}}], "der" -> Sum[ref[p][l]["der"], {p, {"SS", "SR", "RR"}}]|>, {l, 0, lmax}]];
   ref];
refVal[ref_, p_, l_] := ref[p][l]["val"];

(* ---------------- contour and real-axis nodes for the Green's-function quadrature ---------------- *)
nPan = 8; nGc = 16;
contourNodes[np_, ng_] := Module[{pn = panelNodes[Subdivide[aL, aR, np], ng], a}, Table[{N[yC[SetPrecision[pt[[1]], 50]], 50], N[(D[yC[a], a] /. a -> SetPrecision[pt[[1]], 50]) pt[[2]], 50]}, {pt, pn}]];
realNodes = panelNodes[Join[gIn, {4, 17/4, aL}], 20] ~Join~ panelNodes[Join[{aR, 23/4, 6}, gOut], 20];          (* breakpoints at the exterior grid points; finer near the tube, where the source grows like 1/(r - r0)^2 *)
(* exterior field: -(1/(2l+1)) Re[ v(r) Sum_{Re s < r} w s^2 u(s) S(s) + u(r) Sum_{Re s > r} w s^2 v(s) S(s) ], over real and contour nodes *)
extField[l_, rN_, rS_, cN_, cS_, grid_] := Module[{s = rN[[All, 1]], w = rN[[All, 2]], z = cN[[All, 1]], wz = cN[[All, 2]], us, vs, uz, vz},
   us = uFun[l, #] & /@ s; vs = vFun[l, #] & /@ s; uz = uFun[l, #] & /@ z; vz = vFun[l, #] & /@ z;
   Table[-(1/(2 l + 1)) Re[vFun[l, r] (Total[w s^2 us rS UnitStep[r - s]] + If[r > r0, Total[wz z^2 uz cS], 0]) +
       uFun[l, r] (Total[w s^2 vs rS (1 - UnitStep[r - s])] + If[r < r0, Total[wz z^2 vz cS], 0])], {r, grid}]];
scaledErr[num_, ref_] := With[{mx = Max[Abs[ref]]}, If[mx == 0, Max[Abs[num - ref]], Max[Abs[num - ref]]/mx]];      (* absolute if the reference vanishes (RR for l >= 5) *)

(* ---------------- closed-form continued sources ---------------- *)
blendPair[{in_, out_}] := logOut[in] + (logOut[out] - logOut[in]) HExt[zVar];
cfSS[l_] := cfSS[l] = (1/2) dzzExpr[Function[j, If[j < 0, 0, blendPair[ssUndiffPair["exact", j]]]], l];
cfSR[l_] := cfSR[l] = (1/2) dzzExpr[Function[j, If[j < 0, 0, blendPair[srUndiffPair["exact", j]]]], l];
pwSS[l_] := pwSS[l] = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, logReal[ssUndiffPair["exact", j][[b]]]]], l], {b, 2}];        (* pointwise, both sides *)
pwSR[l_] := pwSR[l] = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, srUndiffPair["exact", j][[b]]]], l], {b, 2}];
pwRR[l_] := pwRR[l] = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, rrUndiffPair["exact", j][[b]]]], l], {b, 2}];
pwAt[pair_, r_] := evR[pair[[If[r < r0, 1, 2]]], r];
(* sources for NDSolve are tabulated on piecewise Chebyshev interpolants (massless Stage 8 technique): the closed-form expressions are
   evaluated once at the nodes (50 digits) and the interpolant is evaluated inside NDSolve; 24 nodes on 40 contour pieces and on
   16 + 16 real pieces (breaks at r0 where needed) keep the sources to ~1e-25 *)
NAc = 24; aBreaks = Subdivide[aL, aR, 40];
chebPieces[f_, breaks_, nn_] := Module[{cs, g}, cs = Table[chebCoef[f /@ chebX[breaks[[k]], breaks[[k + 1]], nn]], {k, Length[breaks] - 1}];
   g[x_?NumericQ] := Module[{k = Min[Length[cs], Max[1, LengthWhile[breaks, # <= x &]]]}, chebEval[cs[[k]], breaks[[k]], breaks[[k + 1]], x]]; g];      (* NumericQ guard: NDSolve first calls the source symbolically *)
contourCheb[expr_] := Module[{f = Function[a, evR[expr, N[yC[a], 60]]]}, chebPieces[f, aBreaks, NAc]];
realCheb[pair_, breaksIn_, breaksOut_] := Module[{fi = Function[r, evR[pair[[1]], r]], fo = Function[r, evR[pair[[2]], r]], gi, go},
   gi = chebPieces[fi, breaksIn, NAc]; go = chebPieces[fo, breaksOut, NAc]; Module[{g}, g[r_?NumericQ] := If[r < r0, gi[r], go[r]]; g]];
brExt = {Subdivide[1, aL, 16], Subdivide[aR, 10, 16]}; brFull = {Subdivide[1, r0, 16], Subdivide[r0, 10, 16]};
(* 3D regular field as {T, d_r T, d_r^2 T} at fixed mu, for the SR quadrature *)
phiRT := Module[{z, mu, e}, e = phiR3D[z, mu]; Function @@ {{z, mu}, {e, D[e, z], D[e, {z, 2}]}}];
(* point-supported fit with the physical conditions: a point source at x0 responds with v_l(r) outside and u_l(r) inside, so the error field
   on the two exterior regions is fitted as alpha b1 + beta b2 with b1 = [e^{-m rho}/rho]_l (the two-sided physical Coulomb mode: the field of a
   unit charge is -b1/(4 Pi)) and b2 = d_{r0} b1 (dipole), using one point in each region (r = 9/4 and r = 8) *)
fitPts = {9/4, 8};                      (* both must be exterior grid points *)
fitPoint[l_, err_, pts_] := Module[{R0, pr, b1, b2, M},
   pr = Sqrt[4 Pi/(2 l + 1)] {uFun[l, zVar] vFun[l, R0], uFun[l, R0] vFun[l, zVar]};
   b1 = toFunction[pr /. R0 -> r0, zVar]; b2 = toFunction[D[pr, R0] /. R0 -> r0, zVar];
   M = N[{{b1[pts[[1]]][[1]], b2[pts[[1]]][[1]]}, {b1[pts[[2]]][[2]], b2[pts[[2]]][[2]]}}, 40]; LinearSolve[M, err]];

ssRows = {}; nqRows = {}; srRows = {}; solRows = {}; fitRows = {}; noiseRows = {}; recovered = <||>; srcCurves = <||>;
Do[SetCMParameters[<|"m" -> m|>]; clearCMCaches[]; Clear[cfSS, cfSR, pwSS, pwSR, pwRR];
 cfSS[l_] := cfSS[l] = (1/2) dzzExpr[Function[j, If[j < 0, 0, blendPair[ssUndiffPair["exact", j]]]], l];
 cfSR[l_] := cfSR[l] = (1/2) dzzExpr[Function[j, If[j < 0, 0, blendPair[srUndiffPair["exact", j]]]], l];
 pwSS[l_] := pwSS[l] = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, logReal[ssUndiffPair["exact", j][[b]]]]], l], {b, 2}];
 pwSR[l_] := pwSR[l] = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, srUndiffPair["exact", j][[b]]]], l], {b, 2}];
 pwRR[l_] := pwRR[l] = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, rrUndiffPair["exact", j][[b]]]], l], {b, 2}];
 Print["\n--- m = ", mStr[m], " ---"];
 ref = loadRef[m];
 Module[{cN = contourNodes[nPan, nGc], zs, h3 = ssH3["exact"], h3SR, cfS, cfR, rS, rR, quadS, quadR, plainS, scaleS, scaleR, extRef, nqList},
  zs = cN[[All, 1]];
  (* real-axis pointwise sources at the real nodes (outside the tube) *)
  rS = Table[pwAt[pwSS[l], #] & /@ realNodes[[All, 1]], {l, 0, lmax}]; rR = Table[pwAt[pwSR[l], #] & /@ realNodes[[All, 1]], {l, 0, lmax}];
  cfS = Table[evR[cfSS[l], #] & /@ zs, {l, 0, lmax}]; cfR = Table[evR[cfSR[l], #] & /@ zs, {l, 0, lmax}];
  Print["    closed-form sources at ", Length[zs], " contour nodes [", tick[], " s]"];
  (* (A) SS by angular quadrature, blended, all l from one set of evaluations per node *)
  angularNodes[80];
  quadS = Module[{B = Table[blendedModes[z, h3, lmax + 2], {z, zs}]}, Table[(1/2) dzzArr[B[[k]], zs[[k]]][[l + 1, 1]], {l, 0, lmax}, {k, Length[zs]}]];
  plainS = Module[{B = Table[branchModes[z, "out", h3, lmax + 2], {z, zs}]}, Table[(1/2) dzzArr[B[[k]], zs[[k]]][[l + 1, 1]], {l, 0, lmax}, {k, Length[zs]}]];
  Print["    SS quadrature (80 nodes per segment) [", tick[], " s]"];
  Do[Module[{pw = Max[Abs[quadS[[l + 1]] - cfS[[l + 1]]]]/Max[Abs[cfS[[l + 1]]]], pwPlain = Max[Abs[plainS[[l + 1]] - cfS[[l + 1]]]]/Max[Abs[cfS[[l + 1]]]], eQ, eC, eP, rr = refVal[ref, "SS", l]},
     eQ = scaledErr[extField[l, realNodes, rS[[l + 1]], cN, quadS[[l + 1]], gCom], rr];
     eC = scaledErr[extField[l, realNodes, rS[[l + 1]], cN, cfS[[l + 1]], gCom], rr];
     eP = scaledErr[extField[l, realNodes, rS[[l + 1]], cN, plainS[[l + 1]], gCom], rr];
     AppendTo[ssRows, {m, l, pw, eQ, eC, pwPlain, eP}];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (A) SS source by blended complex-radius quadrature = closed-form continuation at the contour nodes (relative)", pw, 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (A) exterior SS field from the quadrature source (contour Green's quadrature) = reference (scaled)", eQ, 10^-11];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (A) exterior SS field from the closed-form contour source = reference (scaled)", eC, 10^-11];
     If[m != 0, csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (A) negative control: one branch alone (no blend) misses the kink and delta content of the massive SS source (exterior error > 1e-4)", If[eP > 10^-4, 0, 1], 0]]],
   {l, lList}];
  (* node-count study *)
  nqList = If[smoke, {40, 80}, {32, 48, 64, 80, 120}];
  Do[angularNodes[nq];
   Module[{B = Table[blendedModes[z, h3, lmax + 2], {z, zs}], qS},
    qS = Table[(1/2) dzzArr[B[[k]], zs[[k]]][[l + 1, 1]], {l, 0, lmax}, {k, Length[zs]}];
    AppendTo[nqRows, {m, nq, 2 nq, Max[Table[Max[Abs[qS[[l + 1]] - cfS[[l + 1]]]]/Max[Abs[cfS[[l + 1]]]], {l, lList}]],
      Max[Table[scaledErr[extField[l, realNodes, rS[[l + 1]], cN, qS[[l + 1]], gCom], refVal[ref, "SS", l]], {l, lList}]]}]], {nq, nqList}];
  angularNodes[80];
  Print["    node-count study [", tick[], " s]"];
  (* (B) SR by quadrature: blend of the undifferentiated product modes 2 q e^{-m rho} PhiR/rho *)
  h3SR = srH3["exact", phiRT];
  quadR = Module[{B = Table[blendedModes[z, h3SR, lmax + 2], {z, zs}]}, Table[(1/2) dzzArr[B[[k]], zs[[k]]][[l + 1, 1]], {l, 0, lmax}, {k, Length[zs]}]];
  Do[Module[{pw = Max[Abs[quadR[[l + 1]] - cfR[[l + 1]]]]/Max[Abs[cfR[[l + 1]]]], eQ, eC, rr = refVal[ref, "SR", l]},
     eQ = scaledErr[extField[l, realNodes, rR[[l + 1]], cN, quadR[[l + 1]], gCom], rr];
     eC = scaledErr[extField[l, realNodes, rR[[l + 1]], cN, cfR[[l + 1]], gCom], rr];
     AppendTo[srRows, {m, l, pw, eQ, eC}];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (B) SR source by blended quadrature of 2 q e^{-m rho} PhiR/rho = closed form at the nodes (relative)", pw, 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (B) exterior SR field from the quadrature source = reference (scaled)", eQ, 10^-11];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (B) exterior SR field from the closed-form contour source = reference (scaled)", eC, 10^-11]], {l, lList}];
  Print["    SR quadrature [", tick[], " s]"];
  (* (D) noise: 1e-8 relative noise on the SS mode samples; noise-induced exterior error only (the systematic part cancels) *)
  Module[{eps = 10^-8, gradedNodes, gS, xi, errC, errR, nReal = 3, cN4 = contourNodes[4, 16], quadS4},
   gradedNodes = Module[{bp = Table[r0 - (aR - r0) (1/5)^k, {k, 0, 6}], bq = Table[r0 + (aR - r0) (1/5)^k, {k, 0, 6}]},
     Join[panelNodes[Append[bp, r0], 8], panelNodes[Prepend[Reverse[bq], r0], 8]]];
   gS = Table[pwAt[pwSS[l], #] & /@ gradedNodes[[All, 1]], {l, 0, lmax}];
   quadS4 = Module[{B = Table[blendedModes[z, h3, lmax + 2], {z, cN4[[All, 1]]}]}, Table[(1/2) dzzArr[B[[k]], cN4[[k, 1]]][[l + 1, 1]], {l, 0, lmax}, {k, Length[cN4]}]];
   Do[SeedRandom[100 l + 7];
    errC = Table[xi = RandomVariate[NormalDistribution[], Length[cN4]]; Max[Abs[extField[l, {{2, 0}}, {0}, cN4, eps xi quadS4[[l + 1]], gCom]]]/Max[Abs[refVal[ref, "SS", l]]], {k, nReal}];
    errR = Table[xi = RandomVariate[NormalDistribution[], Length[gradedNodes]]; Max[Abs[physResp[l, gradedNodes, eps xi gS[[l + 1]], gCom]]]/Max[Abs[refVal[ref, "SS", l]]], {k, nReal}];
    AppendTo[noiseRows, {m, l, Length[cN4], Max[errC], Length[gradedNodes], Max[errR], Max[errR]/Max[errC]}], {l, lList}];
   csCheck["m=" <> mStr[m] <> ": (D) noise 1e-8 on the SS samples: the real-axis graded grid is > 100x more sensitive than the contour (max over l)", If[Max[noiseRows[[All, 7]]] > 100, 0, 1], 0]];
  Print["    noise test [", tick[], " s]"];
  (* (C) three-region NDSolve solves with the physical conditions; RR on the real line; point-supported fits; negative control *)
  Do[Module[{srcS, srcR, cS, cR, solS, solR, solRR, lRR, eS, eR, eRR, eT, eCT, rrS = refVal[ref, "SS", l], rrR = refVal[ref, "SR", l], rrRR = refVal[ref, "RR", l], rrT = refVal[ref, "total", l],
     abS, abR, outPts = fitPts, posOut, pwNoDelta, eNoDelta, eFixed, gam, posS, extS, extR},
     srcS = realCheb[pwSS[l], brExt[[1]], brExt[[2]]]; srcR = realCheb[pwSR[l], brExt[[1]], brExt[[2]]];
     cS = contourCheb[cfSS[l]]; cR = contourCheb[cfSR[l]];
     solS = threeRegionPhys[l, srcS, cS]; solR = threeRegionPhys[l, srcR, cR];
     Module[{rr = realCheb[pwRR[l], brFull[[1]], brFull[[2]]]}, lRR = realLinePhys[l, rr, rr]];
     extS = solS["ext"] /@ gCom; extR = solR["ext"] /@ gCom;
     eS = scaledErr[extS, rrS]; eR = scaledErr[extR, rrR]; eRR = scaledErr[lRR["f"] /@ gCom, rrRR]; eT = scaledErr[extS + extR + lRR["f"] /@ gCom, rrT];
     eCT = scaledErr[extS + extR, rrS + rrR];                       (* the contour pieces together, scaled by their own size (RR dominates the total at m = 1) *)
     posOut = Flatten[Position[gCom, #] & /@ outPts];
     abS = fitPoint[l, (extS - rrS)[[posOut]], outPts]; abR = fitPoint[l, (extR - rrR)[[posOut]], outPts];
     (* negative control: pointwise real-line SR without the delta; adding c_l G(r, r0) restores it *)
     pwNoDelta = Module[{g = realCheb[pwSR[l], brFull[[1]], brFull[[2]]]}, realLinePhys[l, g, g]];
     gam = (1/2) dzzDeltaCoeff[Function[j, srUndiffPair["exact", j]], l];
     eNoDelta = scaledErr[pwNoDelta["f"] /@ gCom, rrR]; eFixed = scaledErr[(pwNoDelta["f"] /@ gCom) + gam (greenPhys[l, #, r0] & /@ gCom), rrR];
     recovered[{m, l}] = <|"SS" -> extS, "SR" -> extR, "RR" -> lRR["f"] /@ gCom, "total" -> extS + extR + lRR["f"] /@ gCom, "ref" -> rrT, "contourS" -> solS["contour"], "contourR" -> solR["contour"]|>;
     AppendTo[solRows, {m, l, eS, eR, eRR, eT, Abs[abS[[1]]] + Abs[abS[[2]]], Abs[abR[[1]]] + Abs[abR[[2]]], eNoDelta, eFixed, gam, Abs[solS["imAR"][[1]]], Abs[solR["imAR"][[1]]], eCT}];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) SS + SR (the contour pieces) = reference, scaled by max|SS + SR|", eCT, 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) SS three-region contour solve (physical conditions) = reference (scaled)", eS, 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) SR three-region contour solve = reference (scaled)", eR, 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) RR real-line solve = reference (scaled)", eRR, 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) total SS + SR + RR = reference (scaled)", eT, 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) SS has no point-supported term (|alpha| + |beta| of the monopole + dipole fit at r = 9/4, 8)", Abs[abS[[1]]] + Abs[abS[[2]]], 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) SR has no point-supported term", Abs[abR[[1]]] + Abs[abR[[2]]], 10^-12];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (C) negative control: pointwise real-line SR without its delta is wrong (> 1e-2), and c_l G(r, r0) restores it", If[eNoDelta > 10^-2 && eFixed < 10^-12, 0, 1], 0];
     Print["    l=", l, ": SS ", sci[eS], ", SR ", sci[eR], ", RR ", sci[eRR], ", total ", sci[eT], "; fits ", sci[Abs[abS[[1]]] + Abs[abS[[2]]]], ", ", sci[Abs[abR[[1]]] + Abs[abR[[2]]]], "  [", tick[], " s]"]],
   {l, lList}];
  (* positive control of the fit: add a unit point charge delta^3(x - x0) to the SS contour source, mode Y_l(0) deltaExt/r0^2: alpha must be -1/(4 Pi) *)
  Module[{l = 2, cS, srcS, sol, ab, posOut},
   srcS = realCheb[pwSS[l], brExt[[1]], brExt[[2]]];
   cS = contourCheb[cfSS[l] + Y0[l] deltaExt[zVar]/r0^2];
   sol = threeRegionPhys[l, srcS, cS]; posOut = Flatten[Position[gCom, #] & /@ fitPts];
   ab = fitPoint[l, ((sol["ext"] /@ gCom) - refVal[ref, "SS", l])[[posOut]], fitPts];
   AppendTo[fitRows, {m, "positive control: + delta^3 at x0 (l = 2)", N[ab[[1]], 12], N[ab[[2]], 12]}];
   csCheck["m=" <> mStr[m] <> ": (C) positive control: an added unit delta^3 at x0 is detected by the fit, alpha = -1/(4 Pi), beta = 0", Abs[ab[[1]] + 1/(4 Pi)] + Abs[ab[[2]]], 10^-12]];
  (* source curves for the figure (l = 2): real axis and contour *)
  srcCurves[m] = Module[{rr = Join[Subdivide[1, 9/2, 60] // Rest, Subdivide[11/2, 10, 60] // Most], rrT = Join[Subdivide[36/10, 49/10, 60], Subdivide[51/10, 64/10, 60]], aa = Subdivide[aL, aR, 80], l = 2},
     <|"real" -> Table[{r, pwAt[pwSS[l], r], pwAt[pwSR[l], r]}, {r, rrT}], "contour" -> Table[{a, evR[cfSS[l], N[yC[a], 40]], evR[cfSR[l], N[yC[a], 40]]}, {a, aa}],
       "SSmodes" -> Table[{r, Sequence @@ Table[pwAt[pwSS[ll], r], {ll, 0, Min[lmax, 4], 2}]}, {r, Join[Subdivide[1, 49/10, 40], Subdivide[51/10, 10, 40]]}]|>];
  Print["    figure data [", tick[], " s]"]],
 {m, mList}];

(* ================= m = 0 against the massless Stage 9 Test 1 ================= *)
If[MemberQ[mList, 0], Module[{ml = {{0, 1.06*^-19, 4.09*^-20, -0.212694462109}, {1, 3.19*^-19, 3.62*^-20, -0.368397614861}, {2, 1.75*^-18, 3.60*^-19, -0.475599275713}, {3, 2.87*^-17, 1.51*^-17, -0.56273665198},
    {4, 4.22*^-17, 1.58*^-16, -0.638083386326}, {5, 8.85*^-17, 2.31*^-16, -0.705427725801}, {6, 1.74*^-16, 3.34*^-16, -0.76688078914}}},
  (* massless results/stage9_test1.csv (Cauchy data at r = 1): SS_err, SR_err, delta_coeff; the PRL numbers with physical conditions are < 3e-18 *)
  Do[Module[{row = Select[solRows, #[[1]] == 0 && #[[2]] == l &][[1]]},
     csCheck["m=0 l=" <> ToString[l] <> ": SR delta coefficient equals the massless value c_l = " <> ToString[ml[[l + 1, 4]]], Abs[row[[11]] - ml[[l + 1, 4]]], 10^-11];
     csCheck["m=0 l=" <> ToString[l] <> ": exact-puncture SS and SR errors at m = 0 are at the massless level (< 1e-14)", Max[row[[3]], row[[4]]], 10^-14]], {l, lList}]]];

(* ================= tables ================= *)
md1 = StringRiffle[Join[{"| m | l | SS quad vs closed form (pointwise rel.) | exterior SS, quad source | exterior SS, closed-form source | one branch, no blend: pointwise | one branch: exterior |", "|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;;]], " | "] <> " |") & /@ ssRows], "\n"];
md2 = StringRiffle[Join[{"| m | nodes per segment | 3D-source evaluations per contour node (two segments) | SS pointwise rel. error vs closed form (max over l) | exterior SS error (max over l) |", "|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> sci[#[[4]]] <> " | " <> sci[#[[5]]] <> " |") & /@ nqRows], "\n"];
md3 = StringRiffle[Join[{"| m | l | SR quad vs closed form (pointwise rel.) | exterior SR, quad source | exterior SR, closed-form source |", "|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;;]], " | "] <> " |") & /@ srRows], "\n"];
md4 = StringRiffle[Join[{"| m | l | SS | SR | RR | total | SS + SR | SS fit abs(alpha)+abs(beta) | SR fit | pointwise SR, no delta | + c_l G | c_l | Im SS(aR) | Im SR(aR) |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;; 6]], " | "] <> " | " <> sci[#[[14]]] <> " | " <> StringRiffle[sci /@ #[[7 ;; 10]], " | "] <> " | " <> ToString[NumberForm[N[#[[11]], 6], 6]] <> " | " <> sci[#[[12]]] <> " | " <> sci[#[[13]]] <> " |") & /@ solRows], "\n"];
md5 = StringRiffle[Join[{"| m | l | contour nodes | noise-induced error, contour | real graded samples | noise-induced error, real axis | ratio |", "|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> sci[#[[4]]] <> " | " <> ToString[#[[5]]] <> " | " <> sci[#[[6]]] <> " | " <> sci[#[[7]]] <> " |") & /@ noiseRows], "\n"];
md6 = StringRiffle[Join[{"| m | control | alpha | beta |", "|---|---|---|---|"}, ("| " <> mStr[#[[1]]] <> " | " <> #[[2]] <> " | " <> ToString[#[[3]]] <> " | " <> ToString[#[[4]]] <> " |") & /@ fitRows], "\n"];
Print["\n", md1, "\n\n", md2, "\n\n", md3, "\n\n", md4, "\n\n", md5, "\n\n", md6];
Export[FileNameJoin[{resDir, "massive_stage3_summary.md"}], "# Massive Stage 3: contour sources and exact-puncture accuracy\n\nScaled exterior errors = max|num - ref|/max|ref| over the 27 exterior grid points in [1, 3.5] U [6.5, 10], each piece scaled by its own reference (RR dominates the total at m = 1, PhiR^2 ~ e^{2 m r}); physical conditions throughout.\n\n## (A) SS by complex-radius angular quadrature (T2)\n\n" <> md1 <>
   "\n\n## node count (T2)\n\n" <> md2 <> "\n\n## (B) SR by angular quadrature\n\n" <> md3 <> "\n\n## (C) three-region contour solves, exact puncture\n\n" <> md4 <> "\n\n## (D) noise (1e-8 relative on the SS samples)\n\n" <> md5 <> "\n\n## fit controls\n\n" <> md6 <> "\n", "Text"];
Export[FileNameJoin[{resDir, "massive_stage3_ss_quadrature.csv"}], Prepend[{mStr[#[[1]]], #[[2]], Sequence @@ (sci /@ #[[3 ;;]])} & /@ ssRows, {"m", "l", "pointwise_rel", "ext_quad", "ext_cf", "plain_pointwise", "plain_ext"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage3_nodes.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], sci[#[[4]]], sci[#[[5]]]} & /@ nqRows, {"m", "nq", "evals_per_node", "pointwise_rel", "ext_err"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage3_sr_quadrature.csv"}], Prepend[{mStr[#[[1]]], #[[2]], Sequence @@ (sci /@ #[[3 ;;]])} & /@ srRows, {"m", "l", "pointwise_rel", "ext_quad", "ext_cf"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage3_solves.csv"}], Prepend[{mStr[#[[1]]], #[[2]], Sequence @@ (sci /@ #[[3 ;; 10]]), N[#[[11]], 12], sci[#[[12]]], sci[#[[13]]], sci[#[[14]]]} & /@ solRows, {"m", "l", "SS", "SR", "RR", "total", "SS_fit", "SR_fit", "SR_no_delta", "SR_delta_restored", "c_l", "Im_SS_aR", "Im_SR_aR", "SS_plus_SR"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage3_noise.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], sci[#[[4]]], #[[5]], sci[#[[6]]], N[#[[7]], 4]} & /@ noiseRows, {"m", "l", "contour_nodes", "err_contour", "real_samples", "err_real", "ratio"}], "CSV"];

(* ================= figures ================= *)
If[! smoke, Module[{mF = 1, sc, pa, pb, pc, pd, pe, rows, cols = <|0 -> Black, 1/5 -> RGBColor[0.1, 0.35, 0.75], 1 -> RGBColor[0.85, 0.45, 0.1]|>, blue = RGBColor[0.1, 0.35, 0.75], orange = RGBColor[0.85, 0.45, 0.1], red = RGBColor[0.8, 0.1, 0.1]},
   (* Fig. 1: contour; l = 2 sources on the real axis (m = 1, with m = 0 dashed) and on the contour *)
   pa = ParametricPlot[{Re[yC[a]], Im[yC[a]]}, {a, aL, aR}, PlotStyle -> Directive[red, AbsoluteThickness[1.4]], Frame -> True, Axes -> False, AspectRatio -> 0.75,
     FrameLabel -> {"Re z", "Im z"}, BaseStyle -> $fs, ImageSize -> 200, PlotRange -> {{3.9, 6.1}, {-0.15, 1.15}},
     Epilog -> {blue, AbsoluteThickness[1.4], Line[{{3.9, 0}, {aL, 0}}], Line[{{aR, 0}, {6.1, 0}}], Black, PointSize[0.02], Point[{r0, 0}], Inset[Style["r\:2080", $fs], {r0, -0.08}]}];
   sc = srcCurves[mF];
   pb = ListPlot[Join[{Select[sc["real"][[All, {1, 2}]], #[[1]] < r0 &], Select[sc["real"][[All, {1, 2}]], #[[1]] > r0 &], Select[sc["real"][[All, {1, 3}]], #[[1]] < r0 &], Select[sc["real"][[All, {1, 3}]], #[[1]] > r0 &]},
       {Select[srcCurves[0]["real"][[All, {1, 2}]], #[[1]] < r0 &], Select[srcCurves[0]["real"][[All, {1, 2}]], #[[1]] > r0 &], Select[srcCurves[0]["real"][[All, {1, 3}]], #[[1]] < r0 &], Select[srcCurves[0]["real"][[All, {1, 3}]], #[[1]] > r0 &]}],
     Joined -> True, PlotStyle -> {Directive[blue, AbsoluteThickness[1.4]], Directive[blue, AbsoluteThickness[1.4]], Directive[orange, AbsoluteThickness[1.4]], Directive[orange, AbsoluteThickness[1.4]],
       Directive[blue, Dashed, AbsoluteThickness[0.8]], Directive[blue, Dashed, AbsoluteThickness[0.8]], Directive[orange, Dashed, AbsoluteThickness[0.8]], Directive[orange, Dashed, AbsoluteThickness[0.8]]},
     Frame -> True, Axes -> False, FrameLabel -> {"r", "\[ScriptL] = 2 source"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> {{3.6, 6.4}, {-3, 3}},
     Prolog -> {GrayLevel[0.92], Rectangle[{aL, -10}, {aR, 10}]}, Epilog -> {Inset[Style["SS (blue), SR (orange); m = 1 solid, m = 0 dashed", $fs], Scaled[{0.03, 0.96}], {Left, Top}]}];
   pc = ListPlot[{sc["contour"][[All, {1, 2}]] /. {a_, v_} :> {a, Re[v]}, sc["contour"][[All, {1, 2}]] /. {a_, v_} :> {a, Im[v]}, sc["contour"][[All, {1, 3}]] /. {a_, v_} :> {a, Re[v]}, sc["contour"][[All, {1, 3}]] /. {a_, v_} :> {a, Im[v]}},
     Joined -> True, PlotStyle -> {Directive[blue, AbsoluteThickness[1.4]], Directive[blue, Dashed, AbsoluteThickness[1.1]], Directive[orange, AbsoluteThickness[1.4]], Directive[orange, Dashed, AbsoluteThickness[1.1]]},
     Frame -> True, Axes -> False, FrameLabel -> {"a", "\[ScriptL] = 2 source on the contour"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     Epilog -> {Inset[Style["m = 1; solid Re, dashed Im", $fs], Scaled[{0.03, 0.96}], {Left, Top}]}];
   Export[FileNameJoin[{figDir, "fig_contour_source.pdf"}], Grid[{{pa, pb, pc}}, Spacings -> 1, Alignment -> Top]];
   Export[FileNameJoin[{figDir, "fig_contour_source.csv"}], Join[{{"curve", "x", "value_re", "value_im"}},
      ({"SS real axis m=1", N[#[[1]], 6], N[#[[2]], 8], 0} & /@ sc["real"]), ({"SR real axis m=1", N[#[[1]], 6], N[#[[3]], 8], 0} & /@ sc["real"]),
      ({"SS real axis m=0", N[#[[1]], 6], N[#[[2]], 8], 0} & /@ srcCurves[0]["real"]), ({"SR real axis m=0", N[#[[1]], 6], N[#[[3]], 8], 0} & /@ srcCurves[0]["real"]),
      ({"SS contour m=1", N[#[[1]], 6], N[Re[#[[2]]], 8], N[Im[#[[2]]], 8]} & /@ sc["contour"]), ({"SR contour m=1", N[#[[1]], 6], N[Re[#[[3]]], 8], N[Im[#[[3]]], 8]} & /@ sc["contour"])], "CSV"];
   (* Fig. 2: SS source modes l = 0, 2, 4 on the real axis for m = 0, 1/5, 1 *)
   pd = ListPlot[Flatten[Table[{Select[srcCurves[m]["SSmodes"][[All, {1, k + 2}]], #[[1]] < r0 &], Select[srcCurves[m]["SSmodes"][[All, {1, k + 2}]], #[[1]] > r0 &]}, {m, {0, 1/5, 1}}, {k, 0, 2}], 2], Joined -> True,
     PlotStyle -> Flatten[Table[{Directive[cols[m], {Dashing[{}], Dashed, Dotted}[[k + 1]], AbsoluteThickness[1.1]], Directive[cols[m], {Dashing[{}], Dashed, Dotted}[[k + 1]], AbsoluteThickness[1.1]]}, {m, {0, 1/5, 1}}, {k, 0, 2}], 2],
     Frame -> True, Axes -> False, FrameLabel -> {"r", "SS source modes"}, BaseStyle -> $fs, ImageSize -> 300, PlotRange -> {{1, 10}, {-2, 2}}, Prolog -> {GrayLevel[0.92], Rectangle[{aL, -10}, {aR, 10}]},
     Epilog -> {Inset[Style["black m = 0, blue 1/5, orange 1; solid \[ScriptL] = 0, dashed 2, dotted 4", $fs], Scaled[{0.03, 0.96}], {Left, Top}]}];
   Export[FileNameJoin[{figDir, "fig_source_modes.pdf"}], pd];
   Export[FileNameJoin[{figDir, "fig_source_modes.csv"}], Join[{{"m", "r", "SS_l0", "SS_l2", "SS_l4"}}, Flatten[Table[Prepend[N[#, 8], mStr[m]] & /@ srcCurves[m]["SSmodes"], {m, {0, 1/5, 1}}], 1]], "CSV"];
   (* Fig. 3: recovered modes (circles) against the reference (lines), m = 1, l = 0, 2, 4; and scaled errors against r *)
   pe = Show[ListPlot[Flatten[Table[{Transpose[{gIn, recovered[{mF, l}]["ref"][[1 ;; Length[gIn]]]}], Transpose[{gOut, recovered[{mF, l}]["ref"][[Length[gIn] + 1 ;;]]}]}, {l, {0, 2, 4}}], 1], Joined -> True,
      PlotStyle -> Flatten[Table[{Directive[c, AbsoluteThickness[1]], Directive[c, AbsoluteThickness[1]]}, {c, {Black, blue, orange}}], 1]],
     ListPlot[Table[Transpose[{gCom, recovered[{mF, l}]["total"]}], {l, {0, 2, 4}}], PlotStyle -> {Black, blue, orange}, PlotMarkers -> {"\[EmptyCircle]", 7}],
     Frame -> True, Axes -> False, FrameLabel -> {"r", "\[CapitalPhi]\!\(\*SubscriptBox[\(2\), \(\[ScriptL]\)]\)"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All, Prolog -> {GrayLevel[0.92], Rectangle[{aL, -10^3}, {aR, 10^3}]},
     Epilog -> {Inset[Style["m = 1: \[ScriptL] = 0 (black), 2 (blue), 4 (orange)", $fs], Scaled[{0.97, 0.96}], {Right, Top}]}];
   rows = Flatten[Table[Module[{e = Abs[recovered[{mF, l}]["total"] - recovered[{mF, l}]["ref"]]/Max[Abs[recovered[{mF, l}]["ref"]]] + 10^-30}, {Transpose[{gIn, e[[1 ;; Length[gIn]]]}], Transpose[{gOut, e[[Length[gIn] + 1 ;;]]}]}], {l, lList}], 1];
   pa = ListLogPlot[rows, Joined -> True, PlotStyle -> Flatten[Table[{#, #} &[Directive[Hue[0.7 l/lmax, 0.8, 0.7], AbsoluteThickness[1]]], {l, lList}], 1], Frame -> True, Axes -> False,
     FrameLabel -> {"r", "scaled error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> {{1, 10}, All}, Prolog -> {GrayLevel[0.92], Rectangle[{aL, 10^-40}, {aR, 1}]},
     Epilog -> {Inset[Style["\[ScriptL] = 0 ... 6", $fs], Scaled[{0.97, 0.96}], {Right, Top}]}];
   Export[FileNameJoin[{figDir, "fig_recovered_modes.pdf"}], Grid[{{pe, pa}}, Spacings -> 1]];
   Export[FileNameJoin[{figDir, "fig_recovered_modes.csv"}], Join[{{"m", "l", "r", "recovered", "reference", "scaled_error"}},
      Flatten[Table[{mStr[m], l, N[gCom[[k]], 6], N[Re[recovered[{m, l}]["total"][[k]]], 12], N[Re[recovered[{m, l}]["ref"][[k]]], 12], N[Abs[recovered[{m, l}]["total"][[k]] - recovered[{m, l}]["ref"][[k]]]/Max[Abs[recovered[{m, l}]["ref"]]], 4]}, {m, mList}, {l, lList}, {k, Length[gCom]}], 2]], "CSV"];
   (* Fig. 4: node-count study *)
   pb = ListLogPlot[Table[Table[{r[[2]], r[[k]]}, {r, Select[nqRows, #[[1]] == m &]}], {m, mList}, {k, {4, 5}}] // Flatten[#, 1] &, Joined -> True, PlotMarkers -> {"\[FilledCircle]", 5},
     PlotStyle -> Flatten[Table[{Directive[cols[m], AbsoluteThickness[1.2]], Directive[cols[m], Dashed, AbsoluteThickness[1.2]]}, {m, mList}], 1], Frame -> True, Axes -> False,
     FrameLabel -> {"Gauss nodes per segment", "SS error"}, BaseStyle -> $fs, ImageSize -> 260, PlotRange -> All,
     Epilog -> {Inset[Style["solid: pointwise vs closed form; dashed: exterior field", $fs], Scaled[{0.97, 0.96}], {Right, Top}], Inset[Style["black m = 0, blue 1/5, orange 1", $fs], Scaled[{0.97, 0.88}], {Right, Top}]}];
   Export[FileNameJoin[{figDir, "fig_ss_nodes.pdf"}], pb];
   Export[FileNameJoin[{figDir, "fig_ss_nodes.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], N[#[[4]], 6], N[#[[5]], 6]} & /@ nqRows, {"m", "nq", "evals_per_node", "pointwise_rel", "ext_err"}], "CSV"]]];
Print["\nStage 3 wall time: ", tick[], " s"];
csSummary[];
exportChecks[FileNameJoin[{resDir, "massive_stage3_checks.csv"}]];
