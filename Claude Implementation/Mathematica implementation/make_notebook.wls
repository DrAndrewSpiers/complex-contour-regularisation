#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Builds ContourSecondOrder_workflow.nb (unevaluated) from the cell list below and checks that it parses.
   Run:  wolframscript -file make_notebook.wls                                                          *)

dir = DirectoryName[$InputFileName];
out = FileNameJoin[{dir, "ContourSecondOrder_workflow.nb"}];

(* string -> input boxes with comments and line breaks kept (uses the front end, started headlessly) *)
toBoxes[str_String] := Module[{b = Quiet@Check[UsingFrontEnd[MathLink`CallFrontEnd[FrontEnd`UndocumentedTestFEParserPacket[str, False]]], $Failed]},
   If[MatchQ[b, {_BoxData, ___}], First[b], str]];
in[str_String] := Cell[toBoxes[str], "Input"];
txt[str_String] := Cell[str, "Text"];

cells = {
  Cell["Complex-contour S/R split: first and second order (toy model)", "Title"],
  txt["Flat 3D Poisson toy model of the complex-contour regularisation project (A. Spiers). Specification: Claude Implementation/singular-field-consistency.tex. All functions live in ContourSecondOrder.wl in this folder; the .wls scripts are the pass/fail test versions of the sections below. Conventions: Y_l = Y_l0 orthonormal, [f]_l = Sqrt[Pi(2l+1)] Integrate[f P_l(mu), {mu,-1,1}], L_l psi = psi'' + 2 psi'/r - l(l+1) psi/r^2, unit-strength first-order mode source delta(r-r0), kappa_l = -(q/r0^2) Sqrt[4 Pi (2l+1)], H(z) = 1 + (I/Pi) Log[z-r0], deltaExt = I/(Pi (z-r0)); real-axis values are boundary values from above. The complex radial variable is the package symbol zVar."],

  Cell["Setup", "Section"],
  in["SetDirectory[NotebookDirectory[]];\nGet[\"ContourSecondOrder.wl\"]"],
  in["SetCSParameters[<|\"q\" -> 1, \"r0\" -> 5, \"c0\" -> 1, \"c1\" -> 1/10|>]  (* rMin = 1, rMax = 10, tube [9/2, 11/2], amp = 1, WP 30 *)"],
  in["ParametricPlot[{Re[yC[a]], Im[yC[a]]}, {a, tubeL[], tubeR[]}, AxesLabel -> {\"Re z\", \"Im z\"}, PlotLabel -> \"contour y(a)\"]"],

  Cell["Stage 1: the first-order singular mode", "Section"],
  txt["psi_p = DSolveValue[L_l psi == deltaExt] /. C[_] -> 0, exactly as in Singular-and-regular-split-works-real-outside-worldtube.nb. DSolve writes the antiderivatives with Log[r0 - z] and Log[z - r0] (l = 0, 1) or Log[r0 - z] and Log[z] (l >= 2). Log[r0 - z] evaluated at real z > r0 with Mathematica's principal branch is the value from BELOW; canonicalizeLogs rewrites it as Log[z - r0] - I Pi, the same analytic function in the upper half plane."],
  in["Table[{l, dsolveSingularMode[l][zVar]}, {l, 0, 3}] // TableForm"],
  in["Table[{l, singularModeCanonical[l][zVar]}, {l, 0, 2}] // TableForm"],
  txt["D_l(r) = Re psi_p(r+i0) - G_l(r) = alpha r^l + beta r^(-l-1).  Result: alpha_0 = r0, beta_0 = r0^2; D_1 = 0; alpha_l = r0^(1-l)/(2l+1), beta_l = 0 for l >= 2 (the note's 'mirror' choice), not the note's beta_l = r0^(l+2)/(2l+1)."],
  in["TableForm[Table[Prepend[singularModeD[l], l], {l, 0, 8}], TableHeadings -> {None, {\"l\", \"alpha\", \"beta\"}}]"],
  in["(* D_l is the same homogeneous solution on both sides of r0 (exact arithmetic) *)\nTable[Simplify[singularModeDValue[2, r] - singularModeD[2] . {r^2, r^-3}], {r, {1, 2, 4, 6, 8, 10}}]"],
  in["With[{l = 2}, Plot[{Re[N[singularModeCanonical[l][r]]], Re[N[dsolveSingularMode[l][r]]], coulombModeUnit[l, r], Re[N[singularModeCorrected[l][r]]]}, {r, 1, 10},\n  PlotLegends -> {\"Re psi_p(r+i0)\", \"Re psi_p naive (principal branch)\", \"G_l\", \"Re phiS corrected\"}, PlotStyle -> {Blue, {Orange, Dashed}, Black, {Red, Dotted}}]]"],
  in["(* (C1) and (C2) for the corrected mode phiS_l = psi_p - D_l *)\nTable[{l, Max[Table[Abs[N[ComplexExpand[Re[singularModeCorrected[l][rr]]] - coulombModeUnit[l, rr], 30]], {rr, {2, 3, 4, 6, 7, 9}}]],\n   Simplify[radialL[l, singularModeCorrected[l][zVar], zVar] - deltaExt[zVar]]}, {l, 0, 6}] // TableForm"],

  Cell["Stage 2: first-order regular field (S_eff = 0)", "Section"],
  txt["The notebook's three-region solve with the exact corrected singular mode. PhiR = c0 + c1 r cos(theta): only l = 0, 1 are nonzero. imAR = {Im phiR(aR), Im phiR'(aR)}, which the notebook discards with Re[valR]."],
  in["res1 = Association@Table[l -> firstOrderRegular[l], {l, 0, 3}];\nTable[{l, res1[l][\"imAR\"]}, {l, 0, 3}]"],
  in["Plot[{Re[res1[1][\"sol2\"][a]], Im[res1[1][\"sol2\"][a]], Re[N[phiRAnalytic[1][yC[a]]]], Im[N[phiRAnalytic[1][yC[a]]]]}, {a, tubeL[], tubeR[]},\n PlotStyle -> {Blue, Red, {Black, Dashed}, {Gray, Dashed}}, PlotLegends -> {\"Re contour\", \"Im contour\", \"Re continuation\", \"Im continuation\"}]"],
  in["(* generic injected field (notebook test), l = 2, phi(1) = 1, phi'(1) = 0 *)\ngen = firstOrderRegular[2, {1, 0}];\nhom = homogeneousWithData[2, 1, {1, 0}];\n{gen[\"imAR\"], Max[Abs[(gen[\"sol3\"][#] - N[hom[#], 30]) & /@ Subdivide[11/2, 10, 20]]]}"],

  Cell["Stage 3: second order, 1/rho^2 case", "Section"],
  txt["Q = \"Product\" (Q[A,B] = A B).  SS and SR: continued source on the contour inside the tube, real source outside; RR: real line throughout. Boundary conditions from the analytic particular solutions (modes of q^2 Log rho, q c0 rho + q c1 (r0 rho + (r cos theta - r0) rho/2), c0^2 r^2/6 + c0 c1 r^3 cos(theta)/5 + c1^2 r^4 cos(theta)^4/12)."],
  in["Q = \"Product\"; l = 2; rMin = csPar[\"rMin\"]; rMax = csPar[\"rMax\"];\nphiRExt = Association@Table[With[{lp = lp}, lp -> Module[{f}, f[r_?NumericQ] := exteriorValue[res1[lp], r]; f]], {lp, 0, 1}];\nphiRTilde = Association@Table[lp -> res1[lp][\"sol2\"], {lp, 0, 1}];\nphiRLine = Association@Table[lp -> realLineSolve[lp, Function[r, 0], {phiRAnalytic[lp][rMin], phiRAnalytic[lp]'[rMin]}], {lp, 0, 1}];"],
  Cell["SS", "Subsection"],
  in["(* Neumann: [1/rho^2]_l = Sqrt[Pi(2l+1)] Q_l(x)/(r r0); continued with 2 Log[z+r0] - 2 Log[z-r0] *)\nssSourceModeExpr[Q, l, zVar, \"continued\"]"],
  in["solSS = pieceAnalyticSolution[Q, \"SS\", l];\nresSS = threeRegionSolve[l, ssSourceReal[Q, l], ssSourceContour[Q, l], {analyticEval[solSS, rMin], analyticEvalD[solSS, rMin]}];\n{resSS[\"imAR\"], Max[Abs[(exteriorValue[resSS, #] - analyticEval[solSS, #]) & /@ Join[Subdivide[1, 9/2, 20], Subdivide[11/2, 10, 20]]]]}"],
  Cell["SR", "Subsection"],
  in["(* K_{l l'} by Gaunt coefficients; Theta -> H continuation; phiR~ from the Stage 2 contour solve *)\nsolSR = pieceAnalyticSolution[Q, \"SR\", l];\nresSR = threeRegionSolve[l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, phiRTilde, \"ThetaToH\"], {analyticEval[solSR, rMin], analyticEvalD[solSR, rMin]}];\nresSRb = threeRegionSolve[l, srSourceReal[Q, l, phiRExt], srSourceContour[Q, l, phiRTilde, \"CorrectedSingular\"], {analyticEval[solSR, rMin], analyticEvalD[solSR, rMin]}];\n{resSR[\"imAR\"], resSRb[\"imAR\"], Max[Abs[(exteriorValue[resSR, #] - analyticEval[solSR, #]) & /@ Join[Subdivide[1, 9/2, 20], Subdivide[11/2, 10, 20]]]]}"],
  Cell["RR", "Subsection"],
  in["solRR = pieceAnalyticSolution[Q, \"RR\", l];\nlineRR = realLineSolve[l, rrSourceReal[Q, l, phiRLine], {analyticEval[solRR, rMin], analyticEvalD[solRR, rMin]}];\nMax[Abs[(lineRR[#] - analyticEval[solRR, #]) & /@ Subdivide[1, 10, 40]]]"],
  Cell["Two-point boundary conditions", "Subsection"],
  in["bvpSS = threeRegionBVP[l, ssSourceReal[Q, l], ssSourceContour[Q, l], {analyticEval[solSS, rMin], analyticEval[solSS, rMax]}];\nMax[Abs[(exteriorValue[bvpSS, #] - analyticEval[solSS, #]) & /@ Join[Subdivide[1, 9/2, 20], Subdivide[11/2, 10, 20]]]]"],
  in["Show[Plot[analyticEval[solSS, r] + analyticEval[solSR, r] + analyticEval[solRR, r], {r, 1, 10}, PlotStyle -> {Black, Dashed}],\n ListPlot[{#, exteriorValue[resSS, #] + exteriorValue[resSR, #] + lineRR[#]} & /@ Join[Subdivide[1, 9/2, 30], Subdivide[11/2, 10, 40]], PlotStyle -> Red],\n PlotLabel -> \"Phi2 mode l = \" <> ToString[l] <> \": analytic (dashed) vs contour method (dots)\"]"],
  txt["Full tables and 3D reconstruction: run stage3b_second_order_solve.wls (results/stage3_*). The Grad.Grad (1/rho^4) case, Q = \"GradDot\", is implemented with the same functions (regular-field inputs are pairs {phi, dphi/dr}); see stage5a/stage5b."]
  };

nb = Notebook[cells, StyleDefinitions -> "Default.nb"];
Export[out, nb, "NB"];
Print["wrote ", out];
(* verify that it parses: Get, Import as NB, and opening it in a (headless) front end *)
chk = Get[out];
Print["Get: head = ", Head[chk], ", cells = ", Length[chk[[1]]], ", Input cells as BoxData: ",
  Count[chk[[1]], Cell[_BoxData, "Input", ___]], " of ", Count[chk[[1]], Cell[_, "Input", ___]]];
imp = Import[out, {"NB", "Notebook"}];
Print["Import[..., {\"NB\", \"Notebook\"}]: head = ", Head[imp], ", cells = ", Length[imp[[1]]]];
fe = Quiet@Check[UsingFrontEnd[Module[{h = NotebookOpen[out, Visible -> False], c}, c = Length[Cells[h]]; NotebookClose[h]; c]], $Failed];
Print["front end NotebookOpen: cells = ", fe];
