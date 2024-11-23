# MultiplierMeromorphic
Implementation in Magma of an algorithm to calculate multiplier ideals of meromorphic functions.

## Requirements
- Install [Magma](https://magma.maths.usyd.edu.au/magma/)
- Download this repository (you only need `MultiplierMeromorphic.m`)
- Download [SingularitiesDim2/](https://github.com/rogolop/SingularitiesDim2)

## Usage example

Calculate the minimal log-resolution of $\frac{f}{g}$ and the multiplier ideals
```math
\mathcal{J}\left(\left(\frac{f}{g}\right)^\lambda\right)
```
for
```math
\begin{align*}
        &f = (y^2-x^3)^4 + x^8 y^5 ,
        \\ &g = y^2-x^3 ,
    \end{align*}
```
in the range $\lambda\in(0,2)$.

### Magma
```
> AttachSpec("./SingularitiesDim2/IntegralClosureDim2.spec");
> Attach("./MultiplierMeromorphic.m");
> P<x, y> := LocalPolynomialRing(RationalField(), 2, "lglex");
> f := (y^2-x^3)^4 + x^8*y^5;
> g := y^2-x^3;
> Nf, Ng, N, Prox := LogResolutionMeromorphic(f, g);
> Js := MultiplierIdealsMeromorphic(f, g : MinJN:=0, MaxJN:=2, ComputeIdeals:=true);
```

### Output
- `Nf`: 
- `Ng`: 
- `N`: 
- `Prox`: 
- `S`: 

