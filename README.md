<div align="center">
  <img src="https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/docs/src/assets/logo_name_300.png">
</div>

BrkgaMpIpr.jl - Julia version
===============================

[![Build Status](https://travis-ci.org/ceandrade/brkga_mp_ipr_julia.svg?branch=master)](https://travis-ci.org/ceandrade/brkga_mp_ipr_julia)

[![Coverage Status](https://coveralls.io/repos/ceandrade/brkga_mp_ipr_julia/badge.svg?branch=master&service=github)](https://coveralls.io/github/ceandrade/brkga_mp_ipr_julia?branch=master)

[![codecov.io](http://codecov.io/github/ceandrade/brkga_mp_ipr_julia/coverage.svg?branch=master)](http://codecov.io/github/ceandrade/brkga_mp_ipr_julia?branch=master)

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
version.**](https://github.com/ceandrade/brkga_mp_ipr_cpp)
We are also developing a
[Python version](https://github.com/ceandrade/brkga_mp_ipr_python)
which is in its earlier stages.
At this moment, we
have no plans to implement the BRKGA-MP-IPR in other languages such as
Java or C#. But if you want to do so, you are must welcome. But
please, keep the API as close as possible to the C++ API (or Julia API in
case you decide go C), and use the best coding and documentation practices of
your chosen language/framework.

If you are not familiar with how BRKGA works, take a look on
[Standard BRKGA](http://dx.doi.org/10.1007/s10732-010-9143-1) and
[Multi-Parent BRKGA](http://dx.doi.org/xxx).
In the future, we will provide a _Prime on BRKGA-MP_
section.

Installation and tests
--------------------------------------------------------------------------------

| **NOTE:** BrkgaMpIpr was developed using Julia 1.2, but it should work fine
on any Julia >= 1.0. Verions <= 0.6 are not supported.|
| --- |

BrkgaMpIpr can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia-repl
pkg> add BrkgaMpIpr
```

BrkgaMpIpr also provides a thorough unit testing that aims to harden and make
the code ready for production environments. From Pkg REPL, just run

```julia-repl
pkg> test
```

!!! note
    The tests take about 10 minutes, mainly because the permutation path relink.

> :warning: **Warning**:
    It is a hard test to test algorithms that use random signals. In
    BrkgaMpIpr, the tests are carefully designed to ensure repeatability. For
    that, we use the Mersenne Twister
    [[1]](https://en.wikipedia.org/wiki/Mersenne_Twister)
    [[2]](http://dx.doi.org/10.1145/272991.272995) as our standard random
    generator number engine, particularly the [version that comes with
    Julia](https://docs.julialang.org/en/v1/stdlib/Random/index.html#Random.MersenneTwister).
    However, it may happen that such engine has slightly different
    implementations across platforms and, therefore, the tests may fail. The
    current version was tested on 64-bit platforms (Mac OS X, GNU/Linux, and
    Windows 10).

Short usage (TL;DR)
--------------------------------------------------------------------------------

The best way to keep it short is to look in the
[`examples`](https://github.com/ceandrade/brkga_mp_ipr_julia/tree/master/examples/tsp) folder
on [the git repo.](https://github.com/ceandrade/brkga_mp_ipr_julia)
From [`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/main_minimal.jl),
which solves the
[Travelling Salesman Problem (TSP)](https://en.wikipedia.org/wiki/Travelling_salesman_problem).
This is a trimmed copy:

```julia
using BrkgaMpIpr

include("tsp_instance.jl")
include("tsp_decoder.jl")

if length(ARGS) < 4
    println("Usage: julia main_minimal.jl <seed> <config-file> " *
            "<num-generations> <tsp-instance-file>")
    exit(1)
end

seed = parse(Int64, ARGS[1])
configuration_file = ARGS[2]
num_generations = parse(Int64, ARGS[3])
instance_file = ARGS[4]

instance = TSP_Instance(instance_file)

brkga_data, control_params = build_brkga(
    instance, tsp_decode!, MINIMIZE, seed, instance.num_nodes,
    configuration_file
)

initialize!(brkga_data)

evolve!(brkga_data, num_generations)

best_cost = get_best_fitness(brkga_data)
@show best_cost
```

You can identify the following basic steps:

1. Create a data structure inherited from `AbstractInstance` to hold
   your input data. This object is passed to the decoder function (example
   [`tsp_instance.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/tsp_instance.jl));

2. Implement a decoder function. This function translates a chromosome (array
   of numbers in the interval [0,1]) to a solution for your problem. The decoder
   must return the solution value or cost to be used as fitness by BRKGA
   (example [`tsp_decoder.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/tsp_decoder.jl));

3. Load the instance and other relevant data;

4. Use `build_brkga` to create a `BrkgaData` that represents
   the internal state of the BRKGA-MP-IPR algorithm;

5. Use `initialize!` to init the BRKGA state;

6. Call `evolve!` to optimize;

7. Call `get_best_fitness` and/or `get_best_chromosome` to
   retrieve the best solution.

[`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/main_minimal.jl)
provides a very minimal example to understand the necessary steps to use the
BRKGA-MP-IPR framework. However,
[`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/main_complete.jl)
provides a full-featured code, handy for scientific use, such as
experimentation and paper writing. This code allows fine-grained control of
the optimization, shows several features of BRKGA-MP-IPR such as the resets,
chromosome injection, and others. It also logs
all optimization steps, _creating outputs easy to be parsed._ **You should use
this code for serious business and experimentation.**

Tutorial and full documentation
--------------------------------------------------------------------------------

Check out the complete tutorial and API documentation:
https://ceandrade.github.io/brkga_mp_ipr_julia

License and Citing
--------------------------------------------------------------------------------

BRKGA-MP-IPR uses a permissive BSD-like license and it can be used as it
pleases you. And since this framework is also part of an academic effort, we
kindly ask you to remember to cite the originating paper of this work. Indeed,
Clause 4 estipulates that "all publications, softwares, or any other materials
mentioning features or use of this software and/or the data used to test it
must cite explicitly the following article":

> C.E. Andrade. R.F. Toso, J.F. Gonçalves, M.G.C. Resende. The Multi-Parent
> Biased Random-key Genetic Algorithm with Implicit Path Relinking. _European
> Jornal of Operational Research_, volume XX, issue X, pages xx-xx, 2019.
> DOI [to be determined](http://dx.doi.org/xxx)

[Check it out the full license.](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/LICENSE.md)

Contributing
--------------------------------------------------------------------------------

[Contribution guidelines for this project](CONTRIBUTING.md)
