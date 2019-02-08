BrkgaMpIpr tutorial and documentation
================================================================================

BrkgaMpIpr.jl provides a _very easy-to-use_ framework for the
Multi-Parent Biased Random-Key Genetic Algorithm with Implict Path Relink
(**BRKGA-MP-IPR**). Assuming that your have a _decoder_ to your problem,
we can setup, run, and extract the value of the best solution in less than
5 commands (obvisiously, you may need few other lines fo code to do a proper
test).

This Julia version provides a framework fast as C/C++, easy-to-code as Python,
and it is much cheaper (indeed, free) than Matlab. But if you are like me and
also like C++, check out the
[**C++ version.**](https://github.com/ceandrade/brkga_mp_ipr).

If you are not familiar with how BRKGA works, take a look on
[Standard BRKGA](http://dx.doi.org/10.1007/s10732-010-9143-1) and
[Multi-Parent BRKGA](http://dx.doi.org/xxx).
In the future, we will provide a _Prime on BRKGA-MP_
section. If you know what _elite set_, _decoder_, and so means,
we can get to the guts on the [Tutorial](@ref).

```@contents
Pages = ["tutorial.md", "api.md"]
```

License and Citing
----------------------------------------

BrkgaMpIpr.jl uses a permissive BSD-like license. Clause 4 estipulates that
"all publications, softwares, or any other materials mentioning features or
use of this software and/or the data used to test it must cite explicitly
the following article":

> C.E. Andrade. R.F. Toso, J.F. Gonçalves, M.G.C. Resende. The Multi-Parent
> Biased Random-key Genetic Algorithm with Implicit Path Relinking. _European
> Jornal of Operational Research_, volume XX, issue X, pages xx-xx, 2019.
> DOI [to be determined](http://dx.doi.org/xxx)
