<div align="center">
  <img src="https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/docs/src/assets/logo_name_300.png">
</div>

BrkgaMpIpr.jl - Julia version
===============================

<table>
<tr>
  <td>Travis Build Status</td>
  <td>
    <a href="https://travis-ci.org/ceandrade/BrkgaMpIpr.jl">
    <img src="https://travis-ci.org/ceandrade/BrkgaMpIpr.jl.svg?branch=master" alt="Build Status" />
    </a>
  </td>
</tr>
<tr>
  <td>AppVeyor Build Status</td>
  <td>
    <a href="https://ci.appveyor.com/project/ceandrade/brkga-mp-ipr-julia">
    <img src="https://ci.appveyor.com/api/projects/status/nea1697qj5jt9sed?svg=true" alt="Build Status" />
    </a>
  </td>
</tr>
<tr>
  <td>Coverage Status</td>
  <td>
    <a href="https://coveralls.io/github/ceandrade/BrkgaMpIpr.jl?branch=master">
    <img src='https://coveralls.io/repos/github/ceandrade/BrkgaMpIpr.jl/badge.svg?branch=master' alt='Coverage Status' /></a>
    </a>
  </td>
</tr>
<tr>
  <td>codecov.io</td>
  <td>
    <a href="http://codecov.io/github/ceandrade/BrkgaMpIpr.jl?branch=master">
    <img src="http://codecov.io/github/ceandrade/BrkgaMpIpr.jl/coverage.svg?branch=master" alt="codecov.io" />
    </a>
  </td>
</tr>
<tr>
  <td>Documentation</td>
  <td>
    <a href="https://ceandrade.github.io/BrkgaMpIpr.jl">
    <img src="https://img.shields.io/badge/Tutorial-API-blue.svg" alt="Documentation" />
    </a>
  </td>
</tr>
<tr>
  <td>License</td>
  <td>
    <a href="https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/LICENSE.md">
    <img src="https://img.shields.io/badge/license-BSD--like-blue" alt="License" />
    </a>
  </td>
</tr>
</table>

_BrkgaMpIpr.jl_ provides a _very easy-to-use_ framework for the
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
Therefore, _BrkgaMpIpr.jl_ should be suitable to be used in production
environments.

If Julia is not suitable to you, we may find useful the
[**C++ version**](https://github.com/ceandrade/brkga_mp_ipr_cpp)
We are also developing a
[**Python version**](https://github.com/ceandrade/brkga_mp_ipr_python)
which is in its earlier stages.
At this moment, we have no plans to implement the BRKGA-MP-IPR in other
languages such as Java or C#. But if you want to do so, you are must welcome.
But please, keep the API as close as possible to the C++ API (or Julia API in
case you decide go C), and use the best coding and documentation practices of
your chosen language/framework.

- [**C++ version**](https://github.com/ceandrade/brkga_mp_ipr_cpp)
- [**Python version**](https://github.com/ceandrade/brkga_mp_ipr_python)

If you are not familiar with how BRKGA works, take a look on
[Standard BRKGA](http://dx.doi.org/10.1007/s10732-010-9143-1) and
[Multi-Parent BRKGA](http://dx.doi.org/xxx).
In the future, we will provide a _Prime on BRKGA-MP_
section.

:computer: Installation and tests
--------------------------------------------------------------------------------

> :information_source: **NOTE:**
    _BrkgaMpIpr.jl_ was developed using Julia 1.2, but it should work fine
    on any Julia >= 1.0. Versions <= 0.6 are not supported.|

_BrkgaMpIpr.jl_ can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia-repl
pkg> add BrkgaMpIpr
```

Or, just use `Pkg` directly:

```julia-repl
julia> import Pkg; Pkg.add("BrkgaMpIpr")
```

_BrkgaMpIpr.jl_ also provides a thorough unit testing that aims to harden and make
the code ready for production environments. From Pkg REPL, just run

```julia-repl
pkg> test BrkgaMpIpr
```

> :information_source: **NOTE:**
    The tests take about 10 minutes, mainly because the permutation path relink.

Although BrkgaMpIpr should work fine on Julia >= 1.2, some tests can fail. This
issue occurs because BrkgaMpIpr uses the JLD package to save the population and
results. JLD uses the HDF5 package, which produces slightly different binaries
of different Julia versions. Although the tests may fail in those cases,
BrkgaMpIpr is functional for regular usage. In the table below, you can see
the testing fails due to JDL binary incompatibility.

| Julia version   | Windows   | Linux   | Mac OS X   |
| :-------------: | :-------: | :-----: | :--------: |
| 1.2             | ![](https://img.shields.io/badge/-Fail-red)   | ![](https://img.shields.io/badge/-Pass-green)   | ![](https://img.shields.io/badge/-Pass-green) |
| 1.3             | ![](https://img.shields.io/badge/-Pass-green) | ![](https://img.shields.io/badge/-Pass-green)   | ![](https://img.shields.io/badge/-Pass-green) |
| 1.4             | ![](https://img.shields.io/badge/-Pass-green) | ![](https://img.shields.io/badge/-Pass-green)   | ![](https://img.shields.io/badge/-Pass-green) |

> :warning: **Warning**:
    Some timing tests may fail when carried out on virtual machines and
    containers. The reason is that in such environments, the code runs much
    slower than on bare metal, and some control loops take much time to finish
    before the time stop.  Usually, the difference is a few seconds, but it is
    enough to break some tests.

> :warning: **Warning**:
    It is a hard test to test algorithms that use random signals. In
    _BrkgaMpIpr.jl_, the tests are carefully designed to ensure repeatability. For
    that, we use the Mersenne Twister
    [[1]](https://en.wikipedia.org/wiki/Mersenne_Twister)
    [[2]](http://dx.doi.org/10.1145/272991.272995) as our standard random
    generator number engine, particularly the [version that comes with
    Julia](https://docs.julialang.org/en/v1/stdlib/Random/index.html#Random.MersenneTwister).
    However, it may happen that such engine has slightly different
    implementations across platforms and, therefore, the tests may fail. The
    current version was tested on 64-bit platforms (Mac OS X, GNU/Linux, and
    Windows 10).

:zap: Usage - TL;DR
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

:books: Tutorial and full documentation
--------------------------------------------------------------------------------

Check out the complete tutorial and API documentation:
https://ceandrade.github.io/BrkgaMpIpr.jl

:black_nib: License and Citing
--------------------------------------------------------------------------------

BRKGA-MP-IPR uses a permissive BSD-like license and it can be used as it
pleases you. And since this framework is also part of an academic effort, we
kindly ask you to remember to cite the originating paper of this work.
Indeed, Clause 4 estipulates that "all publications, softwares, or any other
materials mentioning features or use of this software (as a whole package or
any parts of it) and/or the data used to test it must cite the following
article explicitly:":

> C.E. Andrade. R.F. Toso, J.F. Gonçalves, M.G.C. Resende. The Multi-Parent
> Biased Random-key Genetic Algorithm with Implicit Path Relinking. _European
> Journal of Operational Research_, To appear, 2019.
> DOI https://doi.org/10.1016/j.ejor.2019.11.037

[Check it out the full license.](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/LICENSE.md)

:pencil2: Contributing
--------------------------------------------------------------------------------

[Contribution guidelines for this project](CONTRIBUTING.md)
