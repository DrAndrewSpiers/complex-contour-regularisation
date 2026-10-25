#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Massive scalar, Stage 4 (T1): truncated punctures PhiP_n = q Sum_{k=-1}^{n} (-m)^(k+1) rho^k/(k+1)! used consistently in SS, SR, RR.
   Residual Phi_res = Phi1 - PhiP_n = PhiR + q (e^{-m rho}/rho - PhiP_n/q): first non-polynomial term rho^(n+2) (odd), so Phi_res is C^(n+1)
   in 3D and its modes are C^(n+3) at r0 (Stage 1).  Expected: per-mode contributions ~ l'^-(n+5), error ~ l'max^-(n+4)
   (massless cq rho^n trick: residual C^(n-1), slopes n+2; this is the m = 0 check).
   Hybrid split (PRD Sec. IV.B, k = 0): SR = 2Q[PhiP_n, T0] + 2Q[PhiP_n, Phi_res - T0]; Part 1 closed form on the contour (carries the SR delta),
   Part 2 formed on the real axis and Cauchy-continued; RR on the real axis.  As in the massless Stage 12b every exterior error is the
   Green's-function response (physical conditions) of its error source, Lemma 1:
     H1 (k = 0): truncated SR + RR pointwise - exact (SR_n + RR_n) pointwise, plus (gP1 + gRR - gEx) delta(r - r0);  H0: without gP1;
     product route: H1 plus the shadow term Re S~(x + i0) - pointwise on the tube (l'max <= 20).
   Also: slopes and per-mode powers; tail acceleration (Upton Lambda with the leading power n + 5); the Cauchy step of the contour solve
   (actual contour solves against the distributional line, n = 3); the shadow growth per coupled mode of the product route (Q = AB type).
   Run:  wolframscript -file massive_stage4_convergence.wls [smoke]                                                                      *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourMassive.wl"}]];
resDir = FileNameJoin[{ParentDirectory[dir], "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
smoke = MemberQ[$ScriptCommandLine, "smoke"];
csResetChecks[]; t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
Off[N::meprec]; Off[General::munfl];
r0 = 5; q = 1; c0 = 1; c1 = 1/10; rMin = 1; rMax = 10; aL = 9/2; aR = 11/2;
cases = If[smoke, {{1, 1}}, {{1/5, 1}, {1/5, 3}, {1/5, 5}, {1, 1}, {1, 3}, {1, 5}, {0, 3}}];      (* {m, n}; m = 0 uses the massless trick puncture q/rho + c q rho^n *)
jMax = 8; Mm = If[smoke, 16, 40]; lB = If[smoke, {0, 2}, Range[0, 6]]; Mprod = Min[20, Mm]; lpRoutes = Select[{4, 8, 12, 16, 20, 30, 40}, # <= Mm &];
gIn = Subdivide[1, 7/2, 12]; gOut = Subdivide[13/2, 10, 14]; gCom = Join[gIn, gOut];
gridPlot = Join[Subdivide[1, 9/2, 14], Subdivide[11/2, 10, 18]]; gAll = Union[gCom, gridPlot];
posCom = Flatten[Position[gAll, #] & /@ gCom]; posOut = Flatten[Position[gAll, #] & /@ gOut];
mStr[m_] := ToString[m, InputForm]; mTag[m_] := StringReplace[mStr[m], "/" -> "_"];
caseTag[{m_, n_}] := "m" <> mTag[m] <> "_n" <> ToString[n];
$fs = {FontFamily -> "Times", FontSize -> 9};
Print["=== Massive Stage 4: hybrid-split convergence; cases {m, n} = ", cases, ", l'max <= ", Mm, ", l = ", lB, " ==="];

(* ---------------- nodes and Green's-function responses ---------------- *)
cNodes[nn_, grid_: gCom] := Module[{bp = Union[Join[{1, 4, 9/2, 11/2, 6, r0, 10}, grid]]}, panelNodes[bp, nn]];
nodesB = cNodes[12, gAll];
GphysV[l_, r_, s_List] := -(s^2/(2 l + 1)) (UnitStep[r - s] (uFun[l, #] & /@ s) vFun[l, r] + (1 - UnitStep[r - s]) uFun[l, r] (vFun[l, #] & /@ s));
(* zero-data-at-aL Green's function on the tube: g(aR, s) and its r-derivative *)
g0of[l_, s_] := (s^2/(2 l + 1)) (uFun[l, aR] vFun[l, s] - uFun[l, s] vFun[l, aR]);
g1of[l_, s_] := (s^2/(2 l + 1)) (uFunD[l, aR] vFun[l, s] - uFun[l, s] vFunD[l, aR]);
(* a contour-vs-line data difference (du, du') at aR is the response of a tube source error; with the physical conditions that error gives
   -(1/(2l+1)) v(r) M1 for r > aR and -(1/(2l+1)) u(r) M2 for r < aL, with (du, du') = (1/(2l+1)) {u(aR) M2 - v(aR) M1, u'(aR) M2 - v'(aR) M1} *)
physFromData[l_, {du_, ddu_}, grid_] := Module[{M1, M2},
   {M2, M1} = LinearSolve[N[(1/(2 l + 1)) {{uFun[l, aR], -vFun[l, aR]}, {uFunD[l, aR], -vFunD[l, aR]}}, 40], {du, ddu}];
   Table[-(1/(2 l + 1)) If[r > r0, vFun[l, r] M1, uFun[l, r] M2], {r, grid}]];
csCheck["physFromData reduces to the massless formulas at m = 0 (l = 3)", Module[{l = 3, du = 1/3, ddu = -2/7, m1, m2, mine},
   SetCMParameters[<|"m" -> 0|>]; m1 = aR^-l ((l + 1) du + aR ddu); m2 = -aR^(l + 1) (l du - aR ddu);
   mine = physFromData[l, {du, ddu}, {2, 8}]; Max[Abs[mine - {-m1 2^l, -m2 8^(-l - 1)}/(2 l + 1)]]], 10^-25];

(* ---------------- derivative-list machinery (massless Stage 9) ---------------- *)
(* the closed-form modes of rho^k for L up to 48 cancel ~100 digits at r + r0 = 10 (the puncture is used without a window), so they are
   evaluated from 200-digit input; the result is kept at 50 digits *)
evR2[e_, r_] := N[e /. zVar -> SetPrecision[r, 200], 50];
pwVal[pair_, r_, br_: Automatic] := evR2[pair[[If[br === Automatic, If[r < r0, 1, 2], If[br === "in", 1, 2]]]], r];
gt[l1_, l2_, l3_] := gt[l1, l2, l3] = gaunt[l1, l2, l3];
LamU[l_, nn_] := 2^nn (2 l + 1)/((2 l + nn + 1) Pochhammer[l + nn/2 + 3/2, -nn]);
tModes0[A0_] := {A0 Sqrt[4 Pi]};                                                            (* T_0 = A0: its only mode *)
p1F[A_, A0_, s_] := Table[2 prodD[Sum[gt[j, 0, L] A[[L + 1]], {L, j, j}], Take[{A0 Sqrt[4 Pi], 0, 0, 0}, Length[A[[1]]]]], {j, 0, jMax}];
nodeFcum[{A_, V_}, Mx_] := Module[{srC, rrSh, tr, z3 = 0 V[[1]]},
   srC = Table[2 prodD[Sum[gt[j, lp, L] A[[L + 1]], {L, Abs[j - lp], j + lp}], V[[lp + 1]]], {j, 0, jMax}, {lp, 0, Mx}];
   tr[j_, l1_, l2_] := With[{g = gt[j, l1, l2]}, If[g == 0, z3, g prodD[V[[l1 + 1]], V[[l2 + 1]]]]];
   rrSh = Table[Sum[tr[j, m, l2], {l2, Max[0, m - j], m}] + Sum[tr[j, l1, m], {l1, Max[0, m - j], m - 1}], {j, 0, jMax}, {m, 0, Mx}];
   {Accumulate /@ srC, Accumulate /@ rrSh}];
nodeF[{A_, V_}, M_, rr_: True] := Module[{srC, rrSh, tr, z3 = 0 V[[1]]},
   srC = Table[2 prodD[Sum[gt[j, lp, L] A[[L + 1]], {L, Abs[j - lp], j + lp}], V[[lp + 1]]], {j, 0, jMax}, {lp, 0, M}];
   rrSh = If[rr, Table[Sum[With[{g = gt[j, l1, l2]}, If[g == 0, z3, g prodD[V[[l1 + 1]], V[[l2 + 1]]]]], {l1, 0, M}, {l2, 0, M}], {j, 0, jMax}], None];
   {Total[srC, {2}], If[rr, rrSh, None]}];
pt[l_, Fl_, s_] := (1/2) dzzArr[Fl, s][[l + 1, 1]];

(* ---------------- Cauchy continuation set-up (massless Stage 8) ---------------- *)
cI = {4, 6}; NC = 60; NA = 24; nG = 80; cs = cauchySetup[{cI[[1]], r0, cI[[2]]}, nG];
aBreaks = Subdivide[aL, aR, 40]; aCheb[nn_] := Flatten[Table[chebX[aBreaks[[k]], aBreaks[[k + 1]], nn], {k, Length[aBreaks] - 1}]];
lineNodes = Join[glNodes[aL, r0, 40], glNodes[r0, aR, 40]];
deltaAt[l_, srcVals_, sLine_] := Module[{nb = Length[aBreaks] - 1, cc, f, u, s = lineNodes[[All, 1]], w = lineNodes[[All, 2]]},
   cc = chebCoef /@ Partition[srcVals, Length[srcVals]/nb];
   f[a_?NumericQ] := Module[{k = Min[nb, Max[1, LengthWhile[aBreaks, # <= a &]]]}, chebEval[cc[[k]], aBreaks[[k]], aBreaks[[k + 1]], a]];
   u = solveContourSegment[l, f, {0, 0}];
   {Re[u[aR]], Re[u'[aR]]} - {Total[w (g0of[l, #] & /@ s) sLine], Total[w (g1of[l, #] & /@ s) sLine]}];

(* ================= main loop ================= *)
eB = <||>; eH0 = <||>; eProd = <||>; dB = <||>; scaleP = <||>; scaleT = <||>; slopes = {}; accRows = {}; cauchyRows = {}; growthRows = {}; gammaRows = {}; profRows = {};
loadRef[m_] := Module[{f = FileNameJoin[{resDir, "massive_reference_m" <> mTag[m] <> ".mx"}], ref, ok = False},
   If[FileExistsQ[f], ref = Import[f]; ok = Quiet[TrueQ[AssociationQ[ref] && KeyExistsQ[ref, "total"] && Max[Keys[ref["total"]]] >= Max[lB] &&
        Abs[robinResidual[0, ref["total"][0]["val"][[1]], ref["total"][0]["der"][[1]], rMin]] < 10^-15]]];
   If[! ok, Print["    (reference for m = ", mStr[m], " missing or invalid: recomputing)"]; ref = <||>; Do[ref[p] = exactReference[Max[lB], gCom, {p}], {p, {"SS", "SR", "RR"}}];
    ref["total"] = Association@Table[l -> <|"val" -> Sum[ref[p][l]["val"], {p, {"SS", "SR", "RR"}}], "der" -> Sum[ref[p][l]["der"], {p, {"SS", "SR", "RR"}}]|>, {l, 0, Max[lB]}]];
   ref];
Do[Module[{m = cs0[[1]], n = cs0[[2]], coeffs, ddP, ddRes, nodeAV, nodeAVc, aC, resC, P, AVin, AVout, cin, cout, a0, ep, exPairs, tubeN, FtR, FtC, lpR, LA = jMax + Mm, pwEx},
   SetCMParameters[<|"m" -> m|>]; clearCMCaches[];
   coeffs = If[m == 0, {{-1, 1}, {n, r0^-(n + 1)}}, punctureCoeffs[n]];
   Print["\n--- m = ", mStr[m], ", n = ", n, ": puncture q Sum p_k rho^k, {k, p_k} = ", coeffs, " ---"];
   (* mode expressions and derivative lists (real pairs and Theta -> H continued) *)
   ddP[L_] := ddP[L] = Table[D[#, {zVar, k}] & /@ punctureModePair[coeffs, L], {k, 0, 3}];           (* [[k+1, br]] *)
   ddRes[lp_] := ddRes[lp] = Table[D[#, {zVar, k}] & /@ residualExactPair[coeffs, lp], {k, 0, 3}];
   aC[L_] := aC[L] = Table[D[punctureModeCont[coeffs, L], {zVar, k}], {k, 0, 2}];
   resC[lp_] := resC[lp] = Table[D[residualExactCont[coeffs, lp], {zVar, k}], {k, 0, 2}];
   nodeAV[s_, br_, nd_: 3] := {Table[Table[pwVal[ddP[L][[k + 1]], s, br], {k, 0, nd - 1}], {L, 0, LA}], Table[Table[pwVal[ddRes[lp][[k + 1]], s, br], {k, 0, nd - 1}], {lp, 0, Mm}]};
   nodeAVc[s_] := {Table[evR2[#, s] & /@ aC[L], {L, 0, jMax + Mprod}], Table[evR2[#, s] & /@ resC[lp], {lp, 0, Mprod}]};
   (* exact (SR_n + RR_n) pointwise source pairs: [Phi1^2]_j - [PhiP_n^2]_j, logs cancel; real forms for pointwise evaluation *)
   exPairs[j_] := exPairs[j] = logReal /@ (ssUndiffPair["exact", j] + srUndiffPair["exact", j] + rrUndiffPair["exact", j] - ssUndiffPair[coeffs, j]);
   ep = Function[j, If[j < 0, {0, 0}, ssUndiffPair["exact", j] + srUndiffPair["exact", j] + rrUndiffPair["exact", j] - ssUndiffPair[coeffs, j]]];
   P = Table[Module[{cum = nodeFcum[nodeAV[s, Automatic], Mm]}, Table[(1/2) dzzArr[cum[[1, All, M + 1]] + cum[[2, All, M + 1]], s][[All, 1]], {M, 0, Mm}]], {s, nodesB[[All, 1]]}];
   Print["    truncated sources at ", Length[nodesB], " nodes [", tick[], " s]"];
   AVin = nodeAV[r0, "in", 4]; AVout = nodeAV[r0, "out", 4];
   cin = nodeFcum[AVin, Mm]; cout = nodeFcum[AVout, Mm];
   a0 = Accumulate[Table[Y0[l] pwVal[residualExactPair[coeffs, l], r0, "in"], {l, 0, Mm}]];      (* T_0 = Phi_res_M(x0) by mode sums *)
   csCheck[caseTag[cs0] <> ": mode sum of Phi_res at x0 converges to PhiR(x0) = " <> ToString[N[phiRAtParticle[]]] <> " (l'max = " <> ToString[Mm] <> ", within 1e-2)", Abs[a0[[-1]] - phiRAtParticle[]], 10^-2];
   If[MemberQ[{3, 5}, n] || smoke,
    tubeN = Join[glNodes[aL, r0, 40], glNodes[r0, aR, 40]];
    FtR = Table[nodeF[nodeAV[s, Automatic, 3][[All]] /. {A_, V_} :> {Take[A, jMax + Mprod + 1], Take[V, Mprod + 1]}, Mprod, False], {s, tubeN[[All, 1]]}];
    FtC = Table[nodeF[nodeAVc[s], Mprod, False], {s, tubeN[[All, 1]]}];
    Print["    tube values for the product route cached [", tick[], " s]"]];
   Do[Module[{exE = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, exPairs[j][[b]]]], l], {b, 2}], exPt, gEx, g0, scale, res, ref},
      exPt = Table[pwVal[exE, s], {s, nodesB[[All, 1]]}];
      gEx = (1/2) dzzDeltaCoeff[ep, l]; g0 = greenPhys[l, #, r0] & /@ gAll;
      ref = loadRef[m];
      (* errors are kept in absolute units; two scales: the contour pieces SS + SR (used throughout) and the total (massless convention, m = 0 check only).
         At m = 1 the RR piece dominates the total (PhiR^2 ~ e^{2 m r}, up to 2.8e3 at r = 10), which would deflate every error by 1e5. *)
      scaleP[{m, n, l}] = Max[Abs[ref["SS"][l]["val"] + ref["SR"][l]["val"]]]; scaleT[{m, n, l}] = Max[Abs[ref["total"][l]["val"]]]; scale = 1;
      csCheck[caseTag[cs0] <> " l=" <> ToString[l] <> ": exact (SR_n + RR_n) delta coefficient = c_l = -4 Pi q PhiR(x0) Y_l(0)/r0^2 (the SS_n delta is carried by the closed-form SS_n)", Abs[gEx + 4 Pi q phiRAtParticle[] Y0[l]/r0^2], 10^-25];
      res = Table[Module[{gRR = (1/2) dzzDeltaNum[l, cin[[2, # + 1, M + 1, 2]] &, cout[[2, # + 1, M + 1, 2]] &], gP1, gSR, base},
          gP1 = (1/2) dzzDeltaNum[l, p1F[AVin[[1]], a0[[M + 1]], r0][[# + 1, 2]] &, p1F[AVout[[1]], a0[[M + 1]], r0][[# + 1, 2]] &];
          gSR = (1/2) dzzDeltaNum[l, cin[[1, # + 1, M + 1, 2]] &, cout[[1, # + 1, M + 1, 2]] &];
          If[M == Mm, AppendTo[gammaRows, {m, n, l, gEx, gSR, gP1, gRR}]];
          base = physResp[l, nodesB, P[[All, M + 1, l + 1]] - exPt, gAll];
          {base + (gP1 + gRR - gEx) g0, base + (gRR - gEx) g0, base + (gSR + gRR - gEx) g0}], {M, 0, Mm}];
      eB[{m, n, l}] = res[[All, 1]]; eH0[{m, n, l}] = res[[All, 2]];
      csCheck[caseTag[cs0] <> " l=" <> ToString[l] <> ": H1 (k = 0) = mode truncation (c): Part 1 carries the whole SR delta (max relative difference over l'max)", Max[Abs[res[[All, 1]] - res[[All, 3]]]]/Max[Abs[res[[All, 3]]]], 10^-20];
      csCheck[caseTag[cs0] <> " l=" <> ToString[l] <> ": truncated RR carries no delta (relative to c_l)", Max[Abs[Cases[gammaRows, {m, n, l, _, _, _, g_} :> g]]]/Abs[gEx], 10^-20];
      If[MemberQ[{3, 5}, n] || smoke, Do[Module[{s = tubeN[[All, 1]], w = tubeN[[All, 2]], fb},
          fb = Table[Re[pt[l, FtC[[i]][[1]], s[[i]]]] - pt[l, FtR[[i]][[1]], s[[i]]], {i, Length[s]}];
          eProd[{m, n, l, M}] = (res[[M + 1, 3]] + Table[Total[w GphysV[l, r, s] fb], {r, gAll}])], {M, Select[lpRoutes, # <= Mprod &]}]];
      dB[{m, n, l}] = Prepend[Differences[eB[{m, n, l}]], eB[{m, n, l}][[1]] - (physResp[l, nodesB, -exPt, gAll] - gEx g0)]], {l, lB}];
   Print["    errors done [", tick[], " s]"];
   (* shadow growth per coupled mode (product route, Q = AB type): (b)_l' = response at aR of -2 Im K~P_{l l'} Im phi~res_l' over the tube *)
   If[! smoke && n == 3, Do[Module[{Rv = R, nd, s, w, g0v, g1v, aCv, l = 2, per, growth, Lmax = 30},
       nd = Join[glNodes[r0 - Rv, r0, 30], glNodes[r0, r0 + Rv, 30]]; s = nd[[All, 1]]; w = nd[[All, 2]];
       g0v = (s^2/(2 l + 1)) (uFun[l, r0 + Rv] (vFun[l, #] & /@ s) - (uFun[l, #] & /@ s) vFun[l, r0 + Rv]);
       aCv = Table[evR2[aC[L][[1]], #] & /@ s, {L, 0, l + Lmax}];
       per = Table[Module[{rcv = evR2[resC[lp][[1]], #] & /@ s, kc}, kc = Sum[gt[l, lp, L] aCv[[L + 1]], {L, Abs[l - lp], l + lp}]; Abs[Total[w g0v (-2 Im[kc] Im[rcv])]]], {lp, 0, Lmax}];
       growth = (per[[Lmax + 1]]/per[[Lmax - 4]])^(1/5);
       AppendTo[growthRows, {m, n, Rv, N[growth, 4], N[1 + Rv/r0, 4], per[[{1, 5, 9, 13, 21, 31}]]}]], {R, {1/2, 1/4}}]];
   Print["    growth done [", tick[], " s]"];
   (* Cauchy step: contour solves against the distributional line, H1 k = 0 and H0, n = 3 *)
   If[(n == 3 && m != 0) || smoke, Do[Module[{M = Mv, l = lv, LA2 = jMax + Mv, tCh, aN, yAc, srAt, sCh, sG, sL, FIn, FOut, gSR, lineAR, p1At, p1Ch, p1G, src, gP1, d, pred, du, eph},
       srAt[s_] := pt[l, nodeF[nodeAV[s, Automatic] /. {A_, V_} :> {Take[A, LA2 + 1], Take[V, M + 1]}, M, False][[1]], s];
       FIn = nodeF[AVin /. {A_, V_} :> {Take[A, LA2 + 1], Take[V, M + 1]}, M, False][[1]]; FOut = nodeF[AVout /. {A_, V_} :> {Take[A, LA2 + 1], Take[V, M + 1]}, M, False][[1]];
       gSR = (1/2) dzzDeltaNum[l, FIn[[# + 1, 2]] &, FOut[[# + 1, 2]] &]; lineAR = gSR {g0of[l, r0], g1of[l, r0]};
       tCh = {chebX[cI[[1]], r0, NC], chebX[r0, cI[[2]], NC]}; aN = aCheb[NA]; yAc = N[yC[#], 50] & /@ aN;
       sCh = Map[srAt, tCh, {2}]; sG = Map[srAt, cs["gauss"][[All, All, 1]], {2}]; sL = srAt /@ lineNodes[[All, 1]];
       If[route === "H0", src = cauchyAt[cs, sG, derivCoefs /@ sCh, #] & /@ yAc; gP1 = 0,
        p1At[s_] := pt[l, p1F[nodeAV[s, Automatic][[1]], a0[[M + 1]], s], s];
        p1Ch = Map[p1At, tCh, {2}]; p1G = Map[p1At, cs["gauss"][[All, All, 1]], {2}];
        src = (evR2[(1/2) dzzExpr[Function[j, If[j < 0, 0, 2 Sum[gt[j, 0, L] aC[L][[1]], {L, j, j}] a0[[M + 1]] Sqrt[4 Pi]]], l], #] & /@ yAc) + (cauchyAt[cs, sG - p1G, derivCoefs /@ (sCh - p1Ch), #] & /@ yAc);
        gP1 = (1/2) dzzDeltaNum[l, p1F[AVin[[1]], a0[[M + 1]], r0][[# + 1, 2]] &, p1F[AVout[[1]], a0[[M + 1]], r0][[# + 1, 2]] &]];
       d = deltaAt[l, src, sL] - lineAR; pred = -(gSR - gP1) {g0of[l, r0], g1of[l, r0]};
       du = d - pred; eph = physFromData[l, du, gCom];
       AppendTo[cauchyRows, {m, n, route, M, l, Max[Abs[eph]]/scaleP[{m, n, l}], Max[Abs[eph]]/Max[Abs[If[route === "H0", eH0, eB][{m, n, l}][[M + 1, posCom]]]], Max[Abs[d - pred]]/Max[Abs[lineAR]]}];
       Print["    Cauchy step ", route, " M=", M, " l=", l, ": contour - line - predicted = ", sci /@ (d - pred), " (|line delta| ", sci[Max[Abs[lineAR]]], ")  [", tick[], " s]"]],
      {route, {"H1 k=0", "H0"}}, {Mv, If[smoke, {8}, {8, 20}]}, {lv, If[smoke, {2}, {2, 6}]}]]],
  {cs0, cases}];

(* ================= slopes, acceleration, tables ================= *)
errM[m_, n_, M_] := Max[Table[Max[Abs[eB[{m, n, l}][[M + 1, posCom]]]]/scaleP[{m, n, l}], {l, lB}]];
errMtot[m_, n_, M_] := Max[Table[Max[Abs[eB[{m, n, l}][[M + 1, posCom]]]]/scaleT[{m, n, l}], {l, lB}]];
errAbs[m_, n_, M_] := Max[Table[Max[Abs[eB[{m, n, l}][[M + 1, posCom]]]], {l, lB}]];
errH0[m_, n_, M_] := Max[Table[Max[Abs[eH0[{m, n, l}][[M + 1, posCom]]]]/scaleP[{m, n, l}], {l, lB}]];
errProd[m_, n_, M_] := Max[Table[Max[Abs[eProd[{m, n, l, M}][[posCom]]]]/scaleP[{m, n, l}], {l, lB}]];
expSlope[{m_, n_}] := If[m == 0, n + 2, n + 4];
Do[Module[{m = cs0[[1]], n = cs0[[2]], fitE, lmx, fitD, Mlo = If[smoke, 6, 20]},
   fitE = Fit[Table[{Log[M], Log[errM[m, n, M]]}, {M, Mlo, Mm}], {1, x}, x];
   lmx = lB[[First@Ordering[Table[-Abs[dB[{m, n, l}][[Mm + 1, posOut[[1]]]]], {l, lB}], 1]]];
   fitD = Fit[Table[{Log[lp], Log[Abs[dB[{m, n, lmx}][[lp + 1, posOut[[1]]]]]]}, {lp, If[smoke, 8, 25], Mm}], {1, x}, x];
   AppendTo[slopes, {m, n, -Coefficient[fitE, x], -Coefficient[fitD, x], lmx, errM[m, n, Min[12, Mm]], errM[m, n, Min[20, Mm]], errM[m, n, Mm]}];
   If[! smoke, csCheck[caseTag[cs0] <> ": error slope = " <> ToString[expSlope[cs0]] <> " (n + 4, massive; n + 2 at m = 0) within 0.5 (l'max 20..40)", Abs[-Coefficient[fitE, x] - expSlope[cs0]], 1/2];
    csCheck[caseTag[cs0] <> ": per-mode power = " <> ToString[expSlope[cs0] + 1] <> " within 0.5 (r = 6.5, l' 25..40)", Abs[-Coefficient[fitD, x] - expSlope[cs0] - 1], 1/2]]], {cs0, cases}];
accel[m_, n_, l_, M_, basis_] := Module[{pw = expSlope[{m, n}] + 1, bf, X, Y, coef, tl, wn = If[basis === "Lambda1", 3, 6]},
   bf = Which[basis === "power", Function[lp, {lp^-pw, lp^(-pw - 1)}], basis === "Lambda", Function[lp, {LamU[lp, -pw], LamU[lp, -pw - 2]}], basis === "Lambda1", Function[lp, {LamU[lp, -pw]}]];
   X = N[bf /@ Range[M - wn + 1, M], 30]; Y = dB[{m, n, l}][[M - wn + 2 ;; M + 1]];
   coef = PseudoInverse[X] . Y;
   tl = If[basis === "power", N[{HurwitzZeta[pw, M + 1], HurwitzZeta[pw + 1, M + 1]}, 30], -N[Sum[bf[lp], {lp, 0, M}], 30]];
   eB[{m, n, l}][[M + 1]] + tl . coef];
accM[m_, n_, M_, basis_] := Max[Table[Max[Abs[accel[m, n, l, M, basis][[posCom]]]]/scaleP[{m, n, l}], {l, lB}]];
accRows = Flatten[Table[{cs0[[1]], cs0[[2]], M, errM[cs0[[1]], cs0[[2]], M], accM[cs0[[1]], cs0[[2]], M, "power"], accM[cs0[[1]], cs0[[2]], M, "Lambda"], accM[cs0[[1]], cs0[[2]], M, "Lambda1"]}, {cs0, cases}, {M, 8, Mm}], 1];
(* m = 0 check against the massless Stage 12b (physical conditions): n = 3 slope 4.87, errors 3.90e-5, 3.27e-6, 1.11e-7 at l'max = 12, 20, 40 *)
If[MemberQ[cases, {0, 3}], Module[{s = Select[slopes, #[[1]] == 0 && #[[2]] == 3 &][[1]]},
   csCheck["m=0 n=3 (massless trick puncture): errors at l'max = 12, 20, 40, scaled by max|total| as in the massless test, reproduce the massless Stage 12b (3.90e-5, 3.27e-6, 1.11e-7) within the 3-digit rounding of those numbers (6e-3 relative)", Max[Abs[Table[errMtot[0, 3, M], {M, {12, 20, 40}}]/{3.90*^-5, 3.27*^-6, 1.11*^-7} - 1]], 6 10^-3];
   csCheck["m=0 n=3: slope reproduces the massless 4.87 within 0.02", Abs[s[[3]] - 4.87], 0.02]]];
Do[csCheck[caseTag[{r[[1]], r[[2]]}] <> " " <> r[[3]] <> " M=" <> ToString[r[[4]]] <> " l=" <> ToString[r[[5]]] <> ": Cauchy step of the contour solve (scaled exterior) < 1e-13", r[[6]], 10^-13], {r, cauchyRows}];
Do[csCheck[caseTag[{r[[1]], r[[2]]}] <> " " <> r[[3]] <> " M=" <> ToString[r[[4]]] <> " l=" <> ToString[r[[5]]] <> ": contour - distributional line at aR = predicted delta miss (relative to |c G|) < 1e-8", r[[8]], 10^-8], {r, cauchyRows}];
Print["    slopes and acceleration [", tick[], " s]"];

md1 = StringRiffle[Join[{"| m | n | slope (l'max 20-40) | expected | per-mode power (l' 25-40) | expected | error at l'max = 12 | 20 | 40 | absolute error at 40 | scaled by max abs(total) at 40 |", "|---|---|---|---|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[NumberForm[#[[3]], 3]] <> " | " <> ToString[expSlope[#[[1 ;; 2]]]] <> " | " <> ToString[NumberForm[#[[4]], 3]] <> " | " <> ToString[expSlope[#[[1 ;; 2]]] + 1] <> " | " <> StringRiffle[sci /@ #[[6 ;; 8]], " | "] <> " | " <> sci[errAbs[#[[1]], #[[2]], Mm]] <> " | " <> sci[errMtot[#[[1]], #[[2]], Mm]] <> " |") & /@ slopes], "\n"];
md2 = StringRiffle[Join[{"| m | n | l'max | plain | + power (2 terms) | + Lambda (2 terms) | + Lambda (1 term) | gain Lambda 2 | gain Lambda 1 |", "|---|---|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> StringRiffle[sci /@ #[[4 ;; 7]], " | "] <> " | " <> ToString[NumberForm[#[[4]]/#[[6]], 3]] <> " | " <> ToString[NumberForm[#[[4]]/#[[7]], 3]] <> " |") & /@ Select[accRows, MemberQ[{12, 20, 30, 40}, #[[3]]] &]], "\n"];
routeRows = Flatten[Table[Join[If[M <= Mprod && KeyExistsQ[eProd, {cs0[[1]], cs0[[2]], lB[[1]], M}], {{cs0[[1]], cs0[[2]], "product", M, errProd[cs0[[1]], cs0[[2]], M]}}, {}],
     {{cs0[[1]], cs0[[2]], "H0", M, errH0[cs0[[1]], cs0[[2]], M]}, {cs0[[1]], cs0[[2]], "H1 = (c)", M, errM[cs0[[1]], cs0[[2]], M]}},
     If[M >= 12, {{cs0[[1]], cs0[[2]], "H1 + tail", M, Cases[accRows, {cs0[[1]], cs0[[2]], M, _, _, _, a_} :> a][[1]]}}, {}]], {cs0, cases}, {M, lpRoutes}], 2];
md3 = StringRiffle[Join[{"| m | n | route | " <> StringRiffle[ToString /@ lpRoutes, " | "] <> " |", "|---|---|---|" <> StringRepeat["---|", Length[lpRoutes]]},
    Flatten@Table["| " <> mStr[cs0[[1]]] <> " | " <> ToString[cs0[[2]]] <> " | " <> rt <> " | " <> StringRiffle[Function[Mv, If[# === {}, "-", sci[#[[1]]]] &[Cases[routeRows, {cs0[[1]], cs0[[2]], rt, Mv, e_} :> e]]] /@ lpRoutes, " | "] <> " |",
      {cs0, cases}, {rt, {"product", "H0", "H1 = (c)", "H1 + tail"}}]], "\n"];
md4 = StringRiffle[Join[{"| m | n | route | l'max | l | Cauchy step (scaled exterior) | relative to the route's error | contour - line vs predicted (rel.) |", "|---|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> #[[3]] <> " | " <> ToString[#[[4]]] <> " | " <> ToString[#[[5]]] <> " | " <> sci[#[[6]]] <> " | " <> sci[#[[7]]] <> " | " <> sci[#[[8]]] <> " |") & /@ cauchyRows], "\n"];
md5 = StringRiffle[Join[{"| m | n | R | growth per l' near l' = 30 | 1 + R/r0 | abs(b_l') at l' = 0, 4, 8, 12, 20, 30 |", "|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]], InputForm] <> " | " <> ToString[#[[4]]] <> " | " <> ToString[#[[5]]] <> " | " <> StringRiffle[sci /@ #[[6]], ", "] <> " |") & /@ growthRows], "\n"];
md6 = StringRiffle[Join[{"| m | n | l | c_l exact | truncated SR delta (l'max) | Part 1 delta | truncated RR delta |", "|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> ToString[NumberForm[N[#[[4]], 8], 8]] <> " | " <> ToString[NumberForm[N[#[[5]], 8], 8]] <> " | " <> ToString[NumberForm[N[#[[6]], 8], 8]] <> " | " <> sci[#[[7]]] <> " |") & /@ gammaRows], "\n"];
md7 = StringRiffle[Join[{"| m | n | l'max | max error, r <= 3.5 | r in [3.5, 4.5] | r in [5.5, 6.5] | r in [6.5, 10] |", "|---|---|---|---|---|---|---|"},
    Flatten@Table[Module[{prof = Table[Max[Table[Abs[eB[{cs0[[1]], cs0[[2]], l}][[M + 1, k]]]/scaleP[{cs0[[1]], cs0[[2]], l}], {l, lB}]], {k, Length[gAll]}], regMax},
       regMax[{lo_, hi_}] := Max[Table[If[lo <= gAll[[k]] <= hi, prof[[k]], 0], {k, Length[gAll]}]];
       "| " <> mStr[cs0[[1]]] <> " | " <> ToString[cs0[[2]]] <> " | " <> ToString[M] <> " | " <> StringRiffle[sci[regMax[#]] & /@ {{1, 7/2}, {7/2, 9/2}, {11/2, 13/2}, {13/2, 10}}, " | "] <> " |"], {cs0, cases}, {M, Select[{12, 20}, # <= Mm &]}]], "\n"];
Print["\n", md1, "\n\n", md2, "\n\n", md3, "\n\n", md4, "\n\n", md5, "\n\n", md6, "\n\n", md7];
Export[FileNameJoin[{resDir, "massive_stage4_summary.md"}], "# Massive Stage 4: hybrid-split convergence with truncated punctures (T1)\n\nScaled exterior error = max|error| / max|SS + SR reference| over the 27 exterior grid points, max over l; physical conditions.  (The total is dominated by RR at m = 1, PhiR^2 ~ e^{2 m r}; the massless convention max|total| is used only for the m = 0 comparison.)\n\n## slopes\n\n" <> md1 <>
   "\n\n## tail acceleration (Lambda fits with the leading power n + 5, or n + 3 at m = 0)\n\n" <> md2 <> "\n\n## routes\n\n" <> md3 <> "\n\n## Cauchy step of the contour solve\n\n" <> md4 <>
   "\n\n## shadow growth per coupled mode (product route, Q = AB type, l = 2)\n\n" <> md5 <> "\n\n## delta coefficients\n\n" <> md6 <> "\n\n## error against r\n\n" <> md7 <> "\n", "Text"];
Export[FileNameJoin[{resDir, "massive_stage4_convergence.csv"}], Prepend[Flatten[Table[{mStr[cs0[[1]]], cs0[[2]], M, sci[errM[cs0[[1]], cs0[[2]], M]], sci[errH0[cs0[[1]], cs0[[2]], M]], sci[errAbs[cs0[[1]], cs0[[2]], M]], sci[errMtot[cs0[[1]], cs0[[2]], M]]}, {cs0, cases}, {M, 0, Mm}], 1], {"m", "n", "lpmax", "H1_err", "H0_err", "H1_abs", "H1_scaled_total"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_slopes.csv"}], Prepend[{mStr[#[[1]]], #[[2]], N[#[[3]], 5], N[#[[4]], 5], #[[5]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]]} & /@ slopes, {"m", "n", "slope", "permode_power", "l_of_fit", "err12", "err20", "err40"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_permode.csv"}], Prepend[Flatten[Table[{mStr[cs0[[1]]], cs0[[2]], l, lp, sci[dB[{cs0[[1]], cs0[[2]], l}][[lp + 1, posOut[[1]]]]]}, {cs0, cases}, {l, lB}, {lp, 0, Mm}], 2], {"m", "n", "l", "lp", "d_lp_at_r6.5"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_routes.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], #[[4]], N[#[[5]], 8]} & /@ routeRows, {"m", "n", "route", "lpmax", "max_l_scaled_error"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_acceleration.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]]} & /@ accRows, {"m", "n", "lpmax", "plain", "power_tail", "Lambda_tail", "Lambda1_tail"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_cauchy_step.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], #[[4]], #[[5]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]]} & /@ cauchyRows, {"m", "n", "route", "lpmax", "l", "scaled", "relative_to_route_error", "line_vs_predicted"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_growth.csv"}], Prepend[{mStr[#[[1]]], #[[2]], ToString[#[[3]], InputForm], #[[4]], #[[5]], Sequence @@ (sci /@ #[[6]])} & /@ growthRows, {"m", "n", "R", "growth", "1+R/r0", "b0", "b4", "b8", "b12", "b20", "b30"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_delta.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], N[#[[4]], 12], N[#[[5]], 12], N[#[[6]], 12], sci[#[[7]]]} & /@ gammaRows, {"m", "n", "l", "gamma_exact", "gamma_SR_trunc", "gamma_Part1", "gamma_RR_trunc"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage4_profile.csv"}], Prepend[Flatten[Table[{mStr[cs0[[1]]], cs0[[2]], M, N[gAll[[k]], 6], sci[Max[Table[Abs[eB[{cs0[[1]], cs0[[2]], l}][[M + 1, k]]]/scaleP[{cs0[[1]], cs0[[2]], l}], {l, lB}]]]}, {cs0, cases}, {M, Select[{12, 20, 40}, # <= Mm &]}, {k, Length[gAll]}], 2], {"m", "n", "lpmax", "r", "max_l_err"}], "CSV"];

(* ================= figures ================= *)
If[! smoke, Module[{cols = <|1 -> RGBColor[0.55, 0.3, 0.7], 3 -> RGBColor[0.1, 0.35, 0.75], 5 -> RGBColor[0.15, 0.6, 0.25]|>, pa, pb, pc, pd, fl, casesM, red = RGBColor[0.8, 0.1, 0.1]},
   fl[s_, m_, n_] := Module[{a = Log[errM[m, n, 30]] + s Log[30]}, Table[{M, Exp[a - s Log[M]]}, {M, {15, Mm}}]];
   casesM[m_] := Select[cases, #[[1]] == m &];
   (* (a), (b): error against l'max for m = 1/5 and m = 1 with fitted slopes (m = 0, n = 3 dashed black) *)
   {pa, pb} = Table[ListLogLogPlot[Join[Table[Table[{M, errM[m, n, M]}, {M, 2, Mm}], {n, casesM[m][[All, 2]]}], Table[fl[Select[slopes, #[[1]] == m && #[[2]] == n &][[1, 3]], m, n], {n, casesM[m][[All, 2]]}], {Table[{M, errM[0, 3, M]}, {M, 2, Mm}]}],
      Joined -> True, PlotStyle -> Join[Table[Directive[cols[n], AbsoluteThickness[1.3]], {n, casesM[m][[All, 2]]}], Table[Directive[Black, Dashed, AbsoluteThickness[0.6]], {n, casesM[m][[All, 2]]}], {Directive[Black, Dotted, AbsoluteThickness[1]]}],
      Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]\[Prime]\:2098\:2090\:2093", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
      Epilog -> Join[Table[Inset[Style["n = " <> ToString[s[[2]]] <> ": slope " <> ToString[NumberForm[s[[3]], 3]], $fs, cols[s[[2]]]], Scaled[{0.03, 0.03 + 0.08 (3 - Position[casesM[m][[All, 2]], s[[2]]][[1, 1]])}], {Left, Bottom}], {s, Select[slopes, #[[1]] == m &]}],
        {Inset[Style["m = " <> mStr[m] <> "; dotted: m = 0, n = 3 (massless trick)", $fs], Scaled[{0.97, 0.97}], {Right, Top}]}]], {m, {1/5, 1}}];
   (* (c): per-mode contributions at r = 6.5, m = 1 *)
   pc = ListLogLogPlot[Table[Table[{lp, Abs[dB[{1, s[[2]], s[[5]]}][[lp + 1, posOut[[1]]]]]}, {lp, 2, Mm}], {s, Select[slopes, #[[1]] == 1 &]}], Joined -> True,
     PlotStyle -> Table[Directive[cols[n], AbsoluteThickness[1.3]], {n, casesM[1][[All, 2]]}], Frame -> True, Axes -> False,
     FrameLabel -> {"\[ScriptL]\[Prime]", "|per-mode contribution| at r = 6.5"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     Epilog -> Join[Table[Inset[Style["n = " <> ToString[s[[2]]] <> ": power " <> ToString[NumberForm[s[[4]], 3]], $fs, cols[s[[2]]]], Scaled[{0.03, 0.03 + 0.08 (3 - Position[casesM[1][[All, 2]], s[[2]]][[1, 1]])}], {Left, Bottom}], {s, Select[slopes, #[[1]] == 1 &]}], {Inset[Style["m = 1", $fs], Scaled[{0.97, 0.97}], {Right, Top}]}]];
   (* (d): routes for m = 1, n = 3: product (red), H0 (grey), H1 (blue), H1 + tail (green), m = 1/5 dashed *)
   pd = ListLogPlot[Flatten[Table[{Table[{M, errProd[m, 3, M]}, {M, Select[lpRoutes, # <= Mprod &]}], Table[{M, errH0[m, 3, M]}, {M, lpRoutes}], Table[{M, errM[m, 3, M]}, {M, lpRoutes}], Table[{M, Cases[accRows, {m, 3, M, _, _, _, a_} :> a][[1]]}, {M, Select[lpRoutes, # >= 12 &]}]}, {m, {1, 1/5}}], 1],
     Joined -> True, PlotMarkers -> {"\[FilledCircle]", 5}, PlotStyle -> Flatten[Table[{Directive[red, st, AbsoluteThickness[1.3]], Directive[Gray, st, AbsoluteThickness[1.3]], Directive[RGBColor[0.1, 0.35, 0.75], st, AbsoluteThickness[1.3]], Directive[RGBColor[0.15, 0.6, 0.25], st, AbsoluteThickness[1.3]]}, {st, {Dashing[{}], Dashed}}], 1],
     Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]\[Prime]\:2098\:2090\:2093", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
     Epilog -> {Inset[Style["n = 3; m = 1 solid, m = 1/5 dashed", $fs], Scaled[{0.97, 0.97}], {Right, Top}], Inset[Style["red product, grey H0, blue H1 = (c), green H1 + \[CapitalLambda] tail", $fs], Scaled[{0.97, 0.90}], {Right, Top}]}];
   Export[FileNameJoin[{figDir, "fig_convergence.pdf"}], Grid[{{pa, pb}, {pc, pd}}, Spacings -> {1, 1}]];
   Export[FileNameJoin[{figDir, "fig_convergence.csv"}], Prepend[Flatten[Table[{mStr[cs0[[1]]], cs0[[2]], M, N[errM[cs0[[1]], cs0[[2]], M], 6], N[errH0[cs0[[1]], cs0[[2]], M], 6], If[KeyExistsQ[eProd, {cs0[[1]], cs0[[2]], lB[[1]], M}], N[errProd[cs0[[1]], cs0[[2]], M], 6], "-"]}, {cs0, cases}, {M, 2, Mm}], 1], {"m", "n", "lpmax", "H1", "H0", "product"}], "CSV"]]];
Print["\nStage 4 wall time: ", tick[], " s"];
csSummary[];
exportChecks[FileNameJoin[{resDir, "massive_stage4_checks.csv"}]];
