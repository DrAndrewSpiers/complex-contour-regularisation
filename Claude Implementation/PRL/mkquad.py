# Generates fig_quadrature.tex for Appendix C of prl.tex.
# (a) r plane: the contour and the circle |z| = r0.  (b) mu plane: the singular angle mu*(z) crosses the
# straight path [-1,1]; a viable mu contour must detour around it.  (c) rho plane (schematic, broken axis):
# the pole is fixed at rho = 0; the paths rho_- -> w -> rho_+ avoid it, while the straight path at a_R runs through it.
# Curves use the actual contour (r0 = 5, R = 1/2, eta = 1); (b) and (c) are schematic where the axes are broken.
import sys
r0,R,eta=5,0.5,1.0
def y(a): return complex(a, eta*(1-((a-r0)/R)**2)**6 if abs(a-r0)<R else 0)
mus=lambda z:(z/r0+r0/z)/2
lo,hi=4.6,5.0
for _ in range(80):
    m=(lo+hi)/2
    if abs(y(m))<r0: lo=m
    else: hi=m
ac=(lo+hi)/2
BL,VI,GR,OR,GY="blue!70!black","violet!85!black","green!50!black","orange!90!black","gray!60!black"
samples=[("z_1",4.85,BL),("z_c",ac,VI)]
A=[4.5+k/400 for k in range(401)]
f=lambda p:f"({p[0]:.4f},{p[1]:.4f})"
path=lambda pts:" -- ".join(f(p) for p in pts)
def cross(p,c,s=0.09): return r"\draw[%s, very thick] (%.4f,%.4f) -- (%.4f,%.4f) (%.4f,%.4f) -- (%.4f,%.4f);"%(c,p[0]-s,p[1]-s,p[0]+s,p[1]+s,p[0]-s,p[1]+s,p[0]+s,p[1]-s)
out=[]; W=out.append
lab=r"\node[font=\scriptsize, fill=white, inner sep=1pt, %s] at %s {%s};"
# ================= (a) r plane, isotropic =================
Y1=5.9; sa=3.6; ox,oy=4.32,-0.12
P=lambda z:((complex(z).real-ox)*sa,Y1+(complex(z).imag-oy)*sa)
TY=Y1+5.15
W(r"\node[anchor=north west,font=\small] at (0,%.3f) {(a) complex $r$ plane};"%TY)
W(r"\draw[->] %s -- %s node[above left=-1pt] {$\operatorname{Re} z$};"%(f(P(4.32)),f(P(5.68))))
W(r"\draw[very thick, blue!60!black] %s -- %s %s -- %s;"%(f(P(4.34)),f(P(4.5)),f(P(5.5)),f(P(5.66))))
arc=[complex((25-t*t)**0.5,t) for t in [k*1.18/60 for k in range(61)]]
W(r"\draw[dashed, gray] %s;"%path([P(z) for z in arc]))
W(lab%("gray!40!black, anchor=east, fill=none",f(P(complex(4.855,1.1))),r"$|z|=r_0$"))
W(r"\draw[very thick, red!70!black] %s;"%path([P(y(a)) for a in A]))
W(r"\node[red!70!black, align=left, anchor=west, font=\scriptsize] at %s {contour\\$z=y(a)$};"%f(P(complex(5.16,0.98))))
for x,l in [(4.5,"a_L"),(5.0,"r_0"),(5.5,"a_R")]:
    p=P(x); W(r"\draw (%.4f,%.4f) -- (%.4f,%.4f) node[below] {$%s$};"%(p[0],p[1]+0.08,p[0],p[1]-0.08,l))
W(r"\node[gray!40!black, font=\scriptsize, align=center] at %s {$|z|<r_0$:\\$\operatorname{Im}\mu^*<0$};"%f(P(complex(4.54,0.62))))
W(r"\node[gray!40!black, font=\scriptsize, align=center] at %s {$|z|>r_0$:\\$\operatorname{Im}\mu^*>0$};"%f(P(complex(5.44,0.55))))
pos={"z_1":"right","z_c":"above left","z_2":"right"}
for n,a,c in samples: W(r"\fill[%s] %s circle (2.2pt) node[%s, text=%s] {$%s$};"%(c,f(P(y(a))),pos[n],c,n))
W(r"\fill[%s] %s circle (2.2pt); \fill[%s] %s circle (2.2pt);"%(GY,f(P(4.5)),OR,f(P(5.5))))
# ================= (b) mu plane: -1, broken axis, magnified near 1 =================
xm1=6.55; xb0,xb1=7.5,8.5; mx0=0.975; sbx=230; sby=360; myl=-0.0045; ymax=0.0075
def Q(m):
    m=complex(m); return (xb1+(m.real-mx0)*sbx, Y1+(m.imag-myl)*sby)
ax=Q(0)[1]
X1=Q(1)[0]; Xend=Q(1.0108)[0]
W(r"\node[anchor=north west,font=\small] at (6.15,%.3f) {(b) $\mu=\cos\theta$ plane (schematic; magnified near $\mu=1$)};"%TY)
W(r"\draw (6.25,%.4f) -- (%.4f,%.4f); \draw[dashed] (%.4f,%.4f) -- (%.4f,%.4f); \draw[->] (%.4f,%.4f) -- (%.4f,%.4f) node[above left=-1pt] {$\operatorname{Re}\mu$};"%(ax,xb0,ax,xb0,ax,xb1,ax,xb1,ax,Xend+0.35,ax))
for xx in (xb0,xb1): W(r"\draw (%.4f,%.4f) -- (%.4f,%.4f);"%(xx-0.06,ax-0.12,xx+0.06,ax+0.12))
# straight path
W(r"\draw[line width=2.4pt, black!25] (%.4f,%.4f) -- (%.4f,%.4f); \draw[line width=2.4pt, black!25, dash pattern=on 3pt off 2pt] (%.4f,%.4f) -- (%.4f,%.4f); \draw[line width=2.4pt, black!25] (%.4f,%.4f) -- (%.4f,%.4f);"%(xm1,ax,xb0,ax,xb0,ax,xb1,ax,xb1,ax,X1,ax))
# viable contour at z_c (branch sigma=-1): along the axis, over mu*(z_c), back to the axis, on to mu=1
xc=Q(mus(y(ac)))[0]; hw=0.5
W(r"\draw[%s, thick] (%.4f,%.4f) -- (%.4f,%.4f); \draw[%s, thick, dash pattern=on 3pt off 2pt] (%.4f,%.4f) -- (%.4f,%.4f);"%(VI,xm1,ax,xb0,ax,VI,xb0,ax,xb1,ax))
W(r"\draw[%s, thick] (%.4f,%.4f) -- (%.4f,%.4f) .. controls (%.4f,%.4f) and (%.4f,%.4f) .. (%.4f,%.4f) -- (%.4f,%.4f);"%(VI,xb1,ax,xc-hw,ax,xc-hw,ax+0.8,xc+hw,ax+0.8,xc+hw,ax,X1,ax))
W(r"\node[%s, font=\scriptsize, anchor=south, fill=white, inner sep=1pt] at (%.4f,%.4f) {viable contour at $z_c$};"%(VI,xc,ax+0.66))
W(r"\fill (%.4f,%.4f) circle (1.8pt) node[below=2pt, font=\scriptsize] {$-1$};"%(xm1,ax))
W(r"\fill (%.4f,%.4f) circle (1.8pt) node[below=2pt, font=\scriptsize] {$1$};"%(X1,ax))
W(r"\node[font=\scriptsize, anchor=south west, black!60] at (%.4f,%.4f) {straight path};"%(xm1-0.1,ax+0.06))
# track of mu*
W(r"\draw[dashed, gray, thick] %s;"%path([Q(mus(y(a))) for a in A]))
for a in [4.68,5.3]: W(r"\draw[->, gray, thick] %s -- %s;"%(f(Q(mus(y(a)))),f(Q(mus(y(a+0.02))))))
W(r"\node[gray!40!black, font=\scriptsize, anchor=north] at %s {track of $\mu^*(z)$};"%f(Q(complex(0.9875,-0.0031))))
for a,c,t,an in [(4.5,GY,r"$\mu^*(a_L)$","below right=-1pt"),(5.5,OR,r"$\mu^*(a_R)$","above right=-1pt")]:
    p=Q(mus(y(a))); W(r"\draw[%s, thick, fill=white] %s circle (1.7pt) node[%s, font=\scriptsize, text=%s] {%s};"%(c,f(p),an,c,t))
labb={"z_1":("below",0,-0.08),"z_c":("below",0,-0.08),"z_2":("below",0,-0.08)}
for n,a,c in samples:
    p=Q(mus(y(a))); W(cross(p,c)); L=labb[n]
    W(r"\node[%s, text=%s, font=\scriptsize, inner sep=1pt] at (%.4f,%.4f) {$\mu^*(%s)$};"%(L[0],c,p[0]+L[1],p[1]+L[2],n))
# ================= (c) rho plane: schematic, broken axis =================
Y2=0.0; sx=3.6; sy=2.15; xL0=-1.55; xL1=1.0; brk=0.9; xR0=9.35; sxR=4.6
def Xr(re):
    if re<=xL1+1e-9: return (re-xL0)*sx
    return (xL1-xL0)*sx+brk+(re-xR0)*sxR
def T(r):
    r=complex(r); return (Xr(r.real), Y2+(r.imag+1.18)*sy)
axc=T(0)[1]; xbrk0=Xr(xL1); xbrk1=xbrk0+brk; Xmax=Xr(10.62)
W(r"\node[anchor=north west,font=\small] at (0,%.3f) {(c) $\rho$ plane, branch $\sigma=-1$ (schematic; axis broken between $\rho\approx1$ and $\rho\approx2r_0$)};"%(Y2+5.6))
W(r"\draw (0,%.4f) -- (%.4f,%.4f); \draw[dashed] (%.4f,%.4f) -- (%.4f,%.4f); \draw[->] (%.4f,%.4f) -- (%.4f,%.4f) node[above left=-1pt] {$\operatorname{Re}\rho$};"%(axc,xbrk0,axc,xbrk0,axc,xbrk1,axc,xbrk1,axc,Xmax+0.5,axc))
for xx in (xbrk0,xbrk1): W(r"\draw (%.4f,%.4f) -- (%.4f,%.4f);"%(xx-0.06,axc-0.12,xx+0.06,axc+0.12))
W(r"\draw[->] %s -- %s node[below right=-2pt] {$\operatorname{Im}\rho$};"%(f(T(-1.16j)),f(T(1.18j))))
# physical path at a_L (grey band), then the straight path at a_R (through the pole)
W(r"\draw[line width=3.2pt, black!22] %s -- %s;"%(f(T(0.5)),f(T(9.5))))
W(r"\draw[%s, line width=1.6pt, dotted] %s -- %s;"%(OR,f(T(10.5)),f(T(-0.5))))
# tracks
W(r"\draw[dashed, gray, thick] %s;"%path([T(r0-y(a)) for a in A]))
W(r"\draw[dashed, gray, thick] %s;"%path([T(r0+y(a)) for a in A]))
for a in [4.75,5.2]:
    W(r"\draw[->, gray, thick] %s -- %s;"%(f(T(r0-y(a))),f(T(r0-y(a+0.02)))))
    W(r"\draw[->, gray, thick] %s -- %s;"%(f(T(r0+y(a))),f(T(r0+y(a+0.02)))))
# paths via w
cols=samples+[("a_R",5.5,OR)]
for n,a,c in cols:
    z=y(a); W(r"\draw[%s, thick] %s -- %s -- %s;"%(c,f(T(z+r0)),f(T(-0.5j)),f(T(r0-z))))
for n,a,c in cols:
    z=y(a); W(r"\fill[%s] %s circle (2.1pt); \fill[%s] %s circle (2.1pt);"%(c,f(T(z+r0)),c,f(T(r0-z))))
W(r"\fill[%s] %s circle (2.1pt); \fill[%s] %s circle (2.1pt);"%(GY,f(T(r0-y(4.5))),GY,f(T(r0+y(4.5)))))
W(lab%("anchor=south, text="+GY,f((T(0.5)[0]+0.1,axc+0.1)),r"$\rho_+(a_L)$"))
W(lab%("anchor=north, text="+GY,f((T(9.5)[0],axc-0.12)),r"$\rho_-(a_L)$"))
W(r"\fill %s circle (1.7pt);"%f(T(-0.5j)))
W(r"\node[font=\scriptsize, anchor=north east, inner sep=1pt] at (%.4f,%.4f) {$w$};"%(T(-0.5j)[0]-0.05,T(-0.5j)[1]-0.05))
W(cross(T(0),"black",0.11))
W(r"\node[font=\scriptsize, align=right, anchor=south east, inner sep=1pt] at %s {pole, fixed\\at $\rho=0$};"%f((T(0)[0]-0.12,axc+0.1)))
# labels: rho_-
pm={"z_1":("anchor=south east",-0.08,0.08),"z_c":("anchor=south east",-0.1,0.06),"z_2":("anchor=south west",0.08,0.08),"a_R":("anchor=north",0,-0.12)}
for n,a,c in cols:
    p=T(y(a)+r0); s=pm[n]; W(lab%(s[0]+", text="+c,f((p[0]+s[1],p[1]+s[2])),r"$\rho_-(%s)$"%n))
pp={"z_1":("anchor=north west",0.1,-0.08),"z_c":("anchor=west",0.12,0),"z_2":("anchor=north east",-0.08,-0.06),"a_R":("anchor=south east",-0.06,0.08)}
for n,a,c in cols:
    p=T(r0-y(a)); s=pp[n]; W(lab%(s[0]+", text="+c,f((p[0]+s[1],p[1]+s[2])),r"$\rho_+(%s)$"%n))
W(lab%("anchor=south, text=gray!40!black, fill=none",f((T(10.0+1j)[0],T(10+1j)[1]+0.12)),r"track of $\rho_-=z+r_0$"))
W(lab%("anchor=north west, text=gray!40!black, fill=none",f((T(-1j)[0]+0.15,T(-1j)[1]-0.08)),r"track of $\rho_+=r_0-z$"))
W(r"\node[font=\scriptsize, align=left, anchor=north west, text=%s] at (0.05,%.4f) {straight path at $a_R$\\(dotted) runs\\through $\rho=0$};"%(OR,axc-0.1))
doc=r"""% Appendix C figure (fig:quad). Generated by mkquad.py; compile with pdflatex.
\documentclass[tikz,border=2pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb}
\usetikzlibrary{arrows.meta}
\begin{document}
\begin{tikzpicture}[>=Stealth, font=\footnotesize]
"""+"\n".join(out)+"\n\\end{tikzpicture}\n\\end{document}\n"
open(sys.argv[1],"w").write(doc)
