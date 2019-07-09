BrkgaMpIpr.jl Guide and Documentation
================================================================================

BrkgaMpIpr.jl provides a _very easy-to-use_ framework for the
Multi-Parent Biased Random-Key Genetic Algorithm with Implict Path Relink
(**BRKGA-MP-IPR**). Assuming that your have a _decoder_ to your problem,
we can setup, run, and extract the value of the best solution in less than
5 commands (obvisiously, you may need few other lines fo code to do a proper
test).

This Julia version provides a framework as fast as C/C++, as easy-to-code as
Python, and it is much cheaper (indeed, free) than Matlab. Unit and coverage
tests are fully implemented, and all pseudo-random test data were carefully
crafted to guarantee reproducibility (although it is possible that some tests
fail because of different versions of the random number generator).
Therefore, BrkgaMpIpr.jl should be suitable to be used in production
environments.

If you are like me and also like C++, check out the [**C++
version.**](https://github.com/ceandrade/brkga_mp_ipr) At this moment, we
have no plans to implement the BRKGA-MP-IPR in other languages such as
Python, Java, or C#. But if you want to do so, you are must welcome. But
please, keep the API as close as possible to the C++ API (or Julia API in
case you decide go C), and use the best coding and documentation practices of
your chosen language/framework.

If you are not familiar with how BRKGA works, take a look on
[Standard BRKGA](http://dx.doi.org/10.1007/s10732-010-9143-1) and
[Multi-Parent BRKGA](http://dx.doi.org/xxx).
In the future, we will provide a _Prime on BRKGA-MP_
section. If you know what _elite set_, _decoder_, and so means,
we can get to the guts on the [Guide / Tutorial](@ref).

```@contents
Pages = ["guide.md", "api.md", "contributing.md"]
```

License and Citing
--------------------------------------------------------------------------------

BRKGA-MP-IPR uses a permissive BSD-like license and it can be used as it
pleases you. And since this framework is also part of an academic effort, we
kindly ask you to remember to cite the originating paper of this work. Indeed,
Clause 4 estipulates that "all publications, softwares, or any other materials
mentioning features or use of this software and/or the data used to test it
must cite explicitly the following article":

> C.E. Andrade. R.F. Toso, J.F. Gonçalves, M.G.C. Resende. The Multi-Parent
> Biased Random-key Genetic Algorithm with Implicit Path Relinking and its
> real-world applications.  *European Journal of Operational Research*,
> volume XX, issue X, pages xx-xx, 2019.
> DOI [to be determined](http://dx.doi.org/xxx))

About the logo
--------------------------------------------------------------------------------

The logo is just a play with 3 chromosomes crossing with each other
(multi-parent) during the mating process. The lines also represent solutions
paths that encounter with each other generating new solutions during the
path-relink.
