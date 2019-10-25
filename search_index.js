var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#BrkgaMpIpr.jl-Guide-and-Documentation-1",
    "page": "Home",
    "title": "BrkgaMpIpr.jl Guide and Documentation",
    "category": "section",
    "text": "BrkgaMpIpr.jl provides a very easy-to-use framework for the Multi-Parent Biased Random-Key Genetic Algorithm with Implict Path Relink (BRKGA-MP-IPR). Assuming that your have a decoder to your problem, we can setup, run, and extract the value of the best solution in less than 5 commands (obvisiously, you may need few other lines fo code to do a proper test).This Julia version provides a framework as fast as C/C++, as easy-to-code as Python, and it is much cheaper (indeed, free) than Matlab. Unit and coverage tests are fully implemented, and all pseudo-random test data were carefully crafted to guarantee reproducibility (although it is possible that some tests fail because of different versions of the random number generator). Therefore, BrkgaMpIpr.jl should be suitable to be used in production environments.If you are like me and also like C++, check out the C++ version. At this moment, we have no plans to implement the BRKGA-MP-IPR in other languages such as Python, Java, or C#. But if you want to do so, you are must welcome. But please, keep the API as close as possible to the C++ API (or Julia API in case you decide go C), and use the best coding and documentation practices of your chosen language/framework.If you are not familiar with how BRKGA works, take a look on Standard BRKGA and Multi-Parent BRKGA. In the future, we will provide a Prime on BRKGA-MP section. If you know what elite set, decoder, and so means, we can get to the guts on the Guide / Tutorial.Pages = [\"guide.md\", \"api.md\", \"contributing.md\"]"
},

{
    "location": "#License-and-Citing-1",
    "page": "Home",
    "title": "License and Citing",
    "category": "section",
    "text": "BRKGA-MP-IPR uses a permissive BSD-like license and it can be used as it pleases you. And since this framework is also part of an academic effort, we kindly ask you to remember to cite the originating paper of this work. Indeed, Clause 4 estipulates that \"all publications, softwares, or any other materials mentioning features or use of this software and/or the data used to test it must cite explicitly the following article\":C.E. Andrade. R.F. Toso, J.F. Gonçalves, M.G.C. Resende. The Multi-Parent Biased Random-key Genetic Algorithm with Implicit Path Relinking and its real-world applications.  European Journal of Operational Research, volume XX, issue X, pages xx-xx, 2019. DOI to be determined)"
},

{
    "location": "#About-the-logo-1",
    "page": "Home",
    "title": "About the logo",
    "category": "section",
    "text": "The logo is just a play with 3 chromosomes crossing with each other (multi-parent) during the mating process. The lines also represent solutions paths that encounter with each other generating new solutions during the path-relink."
},

{
    "location": "guide/#",
    "page": "Guide / Tutorial",
    "title": "Guide / Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "guide/#Guide-/-Tutorial-1",
    "page": "Guide / Tutorial",
    "title": "Guide / Tutorial",
    "category": "section",
    "text": ""
},

{
    "location": "guide/#Installation-and-tests-1",
    "page": "Guide / Tutorial",
    "title": "Installation and tests",
    "category": "section",
    "text": "BrkgaMpIpr can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and runpkg> add BrkgaMpIprBrkgaMpIpr also provides a thorough unit testing that aims to harden and make the code ready for production environments. From Pkg REPL, just runpkg> testnote: Note\nThe tests take about 10 minutes, mainly because the permutation path relink.warning: Warning\nIt is a hard test to test algorithms that use random signals. In BrkgaMpIpr, the tests are carefully designed to ensure repeatability. For that, we use the Mersenne Twister [1] [2] as our standard random generator number engine, particularly the version that comes with Julia. However, it may happen that such engine has slightly different implementations across platforms and, therefore, the tests may fail. The current version was tested on 64-bit platforms (Mac OS X, GNU/Linux, and Windows 10)."
},

{
    "location": "guide/#TL;DR-1",
    "page": "Guide / Tutorial",
    "title": "TL;DR",
    "category": "section",
    "text": "The best way to keep it short is to look in the examples folder on the git repo. From main_minimal.jl, you can identify the following basic steps:Create a data structure inherited from AbstractInstance to hold your input data. This object is passed to the decoder function (example tsp_instance.jl);\nImplement a decoder function. This function translates a chromosome (array of numbers in the interval [0,1]) to a solution for your problem. The decoder must return the solution value or cost to be used as fitness by BRKGA (example tsp_decoder.jl);\nLoad the instance and other relevant data;\nUse build_brkga to create a BrkgaData that represents the internal state of the BRKGA-MP-IPR algorithm;\nUse initialize! to init the BRKGA state;\nCall evolve! to optimize;\nCall get_best_fitness and/or get_best_chromosome to retrieve the best solution.These are the basic steps, but I do recommend the reading of this guide."
},

{
    "location": "guide/#Getting-started-1",
    "page": "Guide / Tutorial",
    "title": "Getting started",
    "category": "section",
    "text": "BrkgaMpIpr is pretty simple, and you must provide one required data structure representing the problem instance, and one required decoder function to translate chromosomes to solutions.Before you go further, please take a look at the examples folder in the git repo. We will use parts of that code in this guide. There, we solve the classical Traveling Salesman Problem. Given a set of cities and the distances between them (full weighted undirect graph), one must find a minimum-cost tour among all cities, such that each city is visited only once (i.e., find a Hamiltonian cycle of minimum cost). These are the files:tsp_instance.jl: contains the input data structures and helper functions;\ntsp_decoder.jl: contains the decoder function for TSP;\ngreedy_tour.jl: simple heuristic that computes a greedy tour;\nconfig.conf:\nconfig.conf: example of parameter settings;\nmain_minimal.jl: minimal code useful to understand and test the framework. You should start here! Please take a look on this file before continue this tutorial;\nmain_complete.jl: full-featured code, handy for scientific use, such as experimentation and paper writing. This code allows fine-grained control of the optimization, show several features of BRKGA-MP-IPR such as the path-reliking calls, resets, chromosome injection, and others. It also logs all optimization steps, creating outputs easy to be parsed. You should use this code for serious business and experimentation;\ninstances: folder containing some TSP instances for testing.When you call main_minimal.jl or main_complete.jl without arguments, they show the usage. For example, assuming you are using a terminal:$ julia main_minimal.jl\nUsage: julia main_minimal.jl <seed> <config-file> <num-generations> <tsp-instance-file>\n\n$ julia main_complete.jl\nUsage:\n  main_complete.jl -c <config_file> -s <seed> -r <stop_rule> -a <stop_arg> -t <max_time> -i <instance_file> [--no_evolution]\n  main_complete.jl (-h | --help)note: Note\nmain_complete.jl uses the DocOpt package. Please, install it before run this script.So, this is a possible output whe calling main_minimal.jl:$ julia main_minimal.jl 27000001 config.conf 100 instances/brazil58.dat\nReading data...\nBuilding BRKGA data and initializing...\nEvolving 100 generations...\nbest_cost = 37552.0For main_complete.jl, the output is more verbose, since we want to capture as much information as possible to do some statistical analysis. The output should be something close to this:$ julia main_complete.jl -c config.conf -s 2700001 -r Generations -a 100 -t 60 -i instances/brazil58.dat\n------------------------------------------------------\n> Experiment started at 2019-02-13T18:40:11.789\n> Instance: instances/brazil58.dat\n> Configuration: config.conf\n> Algorithm Parameters:\n>  - population_size 2000\n>  - elite_percentage 0.3\n>  - mutants_percentage 0.15\n>  - num_elite_parents 2\n>  - total_parents 3\n>  - bias_type LOGINVERSE\n>  - num_independent_populations 3\n>  - pr_number_pairs 0\n>  - pr_minimum_distance 0.15\n>  - pr_type PERMUTATION\n>  - pr_selection BESTSOLUTION\n>  - alpha_block_size 1.0\n>  - pr_percentage 1.0\n>  - exchange_interval 200\n>  - num_exchange_indivuduals 2\n>  - reset_interval 600\n> Seed: 2700001\n> Stop rule: GENERATIONS\n> Stop argument: 100\n> Maximum time (s): 60.0\n> Number of parallel threads for decoding: 1\n------------------------------------------------------\n\n[18:40:11.87] Reading TSP data...\nNumber of nodes: 58\n\n[18:40:11.906] Generating initial tour...\nInitial cost: 30774.0\n\n[18:40:11.909] Building BRKGA data...\nNew population size: 580\n\n[18:40:12.092] Initializing BRKGA data...\n\n[18:40:12.247] Warming up...\n\n[18:40:12.771] Evolving...\n* Iteration | Cost | CurrentTime\n* 1 | 30774 | 0.03\n* 34 | 30751 | 0.83\n* 35 | 30507 | 0.85\n* 36 | 30088 | 0.87\n* 38 | 30023 | 0.93\n* 39 | 29882 | 0.95\n* 40 | 29665 | 0.97\n* 41 | 29131 | 1.00\n* 57 | 28221 | 1.38\n* 66 | 28211 | 1.59\n* 83 | 28200 | 2.01\n* 86 | 28129 | 2.08\n* 91 | 28118 | 2.19\n[18:40:15.171] End of optimization\n\nTotal number of iterations: 100\nLast update iteration: 91\nTotal optimization time: 2.40\nLast update time: 2.19\nLarge number of iterations between improvements: 33\nTotal path relink time: 0.00\nTotal path relink calls: 0\nNumber of homogenities: 0\nImprovements in the elite set: 0\nBest individual improvements: 0\n\n% Best tour cost: 28118\n% Best tour: 22 8 1 30 13 40 25 9 32 20 53 50 4 18 44 24 58 5 27 43 12 57 23 54 55 2 41 35 10 52 51 47 49 3 48 39 29 36 17 26 19 6 28 14 37 34 56 46 15 45 33 21 11 16 42 38 31 7\n\nInstance,Seed,NumNodes,TotalIterations,TotalTime,TotalPRTime,PRCalls,NumHomogenities,NumPRImprovElite,NumPrImprovBest,LargeOffset,LastUpdateIteration,LastUpdateTime,Cost\nbrazil58.dat,2700001,58,100,2.40,0.00,0,0,0,0,33,91,2.19,28118I hope by now you got your system set up and running. Let\'s see the essential details on how to use the BrkgaMpIpr."
},

{
    "location": "guide/#First-things-first:-basic-data-structures-and-decoder-function-1",
    "page": "Guide / Tutorial",
    "title": "First things first: basic data structures and decoder function",
    "category": "section",
    "text": "All problem information must be encapsulated in a struct inherited from AbstractInstance. AbstractInstance is an empty-abstract struct required in the signature of the decoder function (described further). For example, assume we want to solve the Traveling Salesman Problem. One possible instance struct could be:struct TSP_Instance <: AbstractInstance\n    num_nodes::Int64\n    distances::Array{Float64}\nendSo, note that we have the number of nodes/citiesnum_nodes, and the distance matrix distances. We may need some additional code to load the instance and to deal with the distance matrix (omitted here, see tsp_instance.jl).The second and most important requirement is the decoder function. The mandatory signature of the decoder isdecode!(chromosome::Array{Float64, 1},\n        problem_instance::AbstractInstance,\n        rewrite::Bool = true)::Float64chromosome is a vector of numbers in the interval [0, 1] to be decoded. problem_instance is the data structure containing information about the problem. Such data is used by the decoder to build a solution. rewrite is an optional argument that indicates if the decoder should rewrite the chromosome, in case of local search / local improvements be performed during the decoder process. This flag is critical if you intend to use the Implicit Path Relink (details on path_relink!). The decoder must return a Float64 that is used as the fitness to rank the chromosomes. In general, fitness is the cost/value of the solution, but you may want to use it to penalize solutions that violate the problem constraints, for example.In our TSP example, we have a very simple decoder that generates a permutation of nodes, and compute the cost of the cycle from that permutation (note the used of function distance that returns the distance between two nodes and it is defined on tsp_instance.jl).function tsp_decode!(chromosome::Array{Float64}, instance::TSP_Instance,\n                     rewrite::Bool = true)::Float64\n\n    permutation = Array{Tuple{Float64, Int64}}(undef, instance.num_nodes)\n    for (index, key) in enumerate(chromosome)\n        permutation[index] = (key, index)\n    end\n\n    sort!(permutation)\n\n    cost = distance(instance, permutation[1][2], permutation[end][2])\n    for i in 1:(instance.num_nodes - 1)\n        cost += distance(instance, permutation[i][2], permutation[i + 1][2])\n    end\n\n    return cost\nendWith the instance data and the decoder ready, we can build the BRKGA data structures and perform the optimization."
},

{
    "location": "guide/#Building-BRKGA-MP-IPR-data-structures-1",
    "page": "Guide / Tutorial",
    "title": "Building BRKGA-MP-IPR data structures",
    "category": "section",
    "text": "BrkgaMpIpr framework revolves over a single data structure called BrkgaData that represents the internal state of the BRKGA-MP-IPR algorithm. Since this structure has no constructor, you must build it using one of the Building functions. There are two build_brkga methods:load the parameters from a file:function build_brkga(\n    problem_instance::AbstractInstance,\n    decode_function!::Function,\n    opt_sense::Sense,\n    seed::Int64,\n    chromosome_size::Int64,\n    configuration_file::String,\n    evolutionary_mechanism_on::Bool = true\n)::Tuple{BrkgaData, ExternalControlParams}load the parameters from a hand-made parameter object:function build_brkga(\n    problem_instance::AbstractInstance,\n    decode_function!::Function,\n    opt_sense::Sense,\n    seed::Int64,\n    chromosome_size::Int64,\n    brkga_params::BrkgaParams,\n    evolutionary_mechanism_on::Bool = true,\n)::BrkgaDataBoth methods require a problem_instance to be used in the decode_function! as commented before.note: Note\nTo date, there is not an easy and clean way to force a function type, as we can do in C and C++, and then, decode_function! is declared as a simple, high-level Function type. However, decode_function! must have the expected signature as explained before.You also must indicate whether you are minimizing or maximizing through optimization Sense.A good seed also must be provided for the (pseudo) random number generator (according to this paper). BrkgaMpIpr uses the Mersenne Twister [1] [2].The chromosome_size also must be given. It indicates the length of each chromosome in the population. In general, this size depends on the instance and how the decoder works.Another common argument is evolutionary_mechanism_on which is enabled by default. When disabled, no evolution is performed. The algorithm only decodes the chromosomes and ranks them. On each generation, all population is replaced excluding the best chromosome. This flag helps on implementations of simple multi-start algorithms.As said before, the difference between the two methods is from where the algorithm\'s hyper-parameters come from. In the first version, the algorithm reads the BRKGA, IPR, and extra parameters from a simple text file that looks like this (see config.conf for detailed example):population_size 2000\nelite_percentage 0.30\nmutants_percentage 0.15\nnum_elite_parents 2\ntotal_parents 3\nbias_type LOGINVERSE\nnum_independent_populations 3\npr_number_pairs 0\npr_minimum_distance 0.15\npr_type PERMUTATION\npr_selection BESTSOLUTION\nalpha_block_size 1.0\npr_percentage 1.0\nexchange_interval 200\nnum_exchange_indivuduals 2\nreset_interval 600When reading such file, the algorithm ignores all blank lines, and lines starting with #. During the building process, the building method creates a BrkgaParams object and a ExternalControlParams object. BrkgaParams contains all hyper-parameters regarding BRKGA and IPR methods and is stored in the brand-new BrkgaData returned by the method. ExternalControlParams are parameters that can be used outside the BRKGA-MP-IPR to control several aspects of the optimization. For instance, interval to apply path relink, reset the population, perform population migration, among others. Although their presence is required on the config file, they are not mandatory to the BRKGA-MP-IPR itself.In the second method, we assume we already have a BrkgaParams object, and we just pass it directly to the function. Note that such param data is deep-copied inside BrkgaData.Let\'s take a look in the example from main_minimal.jl:seed = parse(Int64, ARGS[1])\nconfiguration_file = ARGS[2]\nnum_generations = parse(Int64, ARGS[3])\ninstance_file = ARGS[4]\n\ninstance = TSP_Instance(instance_file)\n\nbrkga_data, control_params = build_brkga(\n    instance, tsp_decode!, MINIMIZE, seed, instance.num_nodes,\n    configuration_file\n)This code gets some arguments from the command line and loads a TSP instance. After that, it builds brkga_data. Note that we pass the instance data and the tsp_decode! function. Since we are looking for cycles of minimum cost, we ask for the algorithm MINIMIZE. The starting seed is also given. Since tsp_decode! considers each chromosome key as a node/city, the length of the chromosome must be the number of nodes, i.e., instance.num_nodes. Finally, we also pass the configuration file.Let\'s take a look in a more elaborated example (main_complete.jl):brkga_params, control_params = load_configuration(configuration_file)\n...\nbrkga_params.population_size = min(brkga_params.population_size,\n                                   10 * instance.num_nodes)\n...\nbrkga_data = build_brkga(instance, tsp_decode!, MINIMIZE, seed,\n                         instance.num_nodes, brkga_params, perform_evolution)Here, we first load the configuration file using the helper function load_configuration. Then, we modify the population size to be the minimum between the original parameter and 10x the number of nodes (making the population size proportional to the instance size is a very common strategy used in genetic algorithms). Then, we call build_brkga using the param data instead of the configuration file (there is an additional parameter perform_evolution to tune the evolution either on or off, not necessary at this moment).Now, we have a BrkgaData which will be used in all other functions during the optimization. Note that we can build several BrkgaData objects using different parameters, decoders, or instance data. These structures can be evolved in parallel and mixed-and-matched at your will. Each one holds a self-contained BRKGA state including populations, fitness information, and a state of the random number generator."
},

{
    "location": "guide/#Initialization-and-Warm-start-solutions-1",
    "page": "Guide / Tutorial",
    "title": "Initialization and Warm-start solutions",
    "category": "section",
    "text": "Before starting the optimization, we need to initialize BrkgaData using initialize! function. This procedure initializes the populations and others data structures of the BRKGA. If an initial population (warm start) is supplied, the initialization method completes the remaining individuals, if they do not exist. This method also performs the initial decoding of the chromosomes. Therefore, depending on the decoder implementation, this can take a while, and you may want to time such procedure. The syntax is pretty straightforward:initialize!(brkga_data)warning: Warning\ninitialize! must be called before any optimization methods.warning: Warning\nBrkgaMpIpr performs the decoding of each chromosome in parallel if multi-thread is enabled. Therefore, we must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using a single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing)."
},

{
    "location": "guide/#Warm-start-solutions-1",
    "page": "Guide / Tutorial",
    "title": "Warm-start solutions",
    "category": "section",
    "text": "One good strategy is to bootstrap the main optimization algorithm with good solutions from fast heuristics [1, 2, 3] or even from relaxations of integer linear programming models [4].To do it, you must set these initial solutions before call initialize!. Since BRKGA-MP-IPR does not know the problem structure, you must encode the warm-start solution as chromosomes (vectors in the interval [0, 1]). In other words, you must do the inverse process that decode! does. For instance, this is a piece of code from main_complete.jl showing this process:initial_cost, initial_tour = greedy_tour(instance)\n...\nkeys = sort(rand(instance.num_nodes))\ninitial_chromosome = zeros(instance.num_nodes)\nfor i in 1:instance.num_nodes\n    initial_chromosome[initial_tour[i]] = keys[i]\nend\n...\nset_initial_population!(brkga_data, [initial_chromosome])\ninitialize!(brkga_data)Here, we create one incumbent solution using the greedy heuristic greedy_tour() found in greedy_tour.jl. It gives us initial_tour which is a sequence of nodes to be visited. In the next four lines, we encode initial_tour. First, we create a vector of sorted random keys. Note that this is the same order that tsp_decode! uses. We then create the initial_chromosome, and fill it up with keys according to the nodes\' order in initial_tour. Finally, we use set_initial_population! to assign the incumbent to the initial population. Note that initial_chromosome in between braces because set_initial_population! takes a vector of chromosomes. See its signature:set_initial_population!(brkga_data::BrkgaData,\n                        chromosomes::Array{Array{Float64, 1}, 1})Indeed, you can have as much warm-start solutions as you like, limited to the size of the population. Just remember:warning: Warning\nset_initial_population! must be called BEFORE initialize!."
},

{
    "location": "guide/#Optimization-time:-evolving-the-population-1",
    "page": "Guide / Tutorial",
    "title": "Optimization time: evolving the population",
    "category": "section",
    "text": "Once all data is set up, it is time to evolve the population and perform other operations like path-relinking, shaking, migration, and others. The call is pretty simple:evolve!(brkga_data::BrkgaData, num_generations::Int64 = 1)evolve! evolves all populations for num_generations.For example, in main_minimal.jl, we just evolve the population for a given number of generations directly and then extract the best solution cost.evolve!(brkga_data, num_generations)\nbest_cost = get_best_fitness(brkga_data)On main_complete.jl, we have fine-grained control on the optimization. There, we have a main loop that evolves the population one generation at a time and performs several operations as to hold the best solution, to check whether it is time for path relink, population reset, among others. The advantage of that code is that we can track all optimization details.warning: Warning\nAgain, the decoding of each chromosome is done in parallel if multi-thread is enabled. Therefore, we must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using a single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing)."
},

{
    "location": "guide/#Accessing-solutions/chromosomes-1",
    "page": "Guide / Tutorial",
    "title": "Accessing solutions/chromosomes",
    "category": "section",
    "text": "Since Julia does not offer encapsulation mechanisms to keep data private within data structures, you can access all chromosomes, fitness, and other data members directly from BrkgaData. However, we do not recommend that, unless you are sure what you are doing. So, BrkgaMpIpr offers some helper functions.Usually, we want to access the best chromosome after some iterations. You can use the companion functions:get_best_fitness(brkga_data::BrkgaData)::Float64get_best_chromosome(brkga_data::BrkgaData)::Array{Float64, 1}get_best_fitness returns the value/fitness of the best chromosome across all populations.get_best_chromosome returns a copy of the best chromosome across all populations. You may want to extract an actual solution from such chromosome, i.e., to apply a decoding function that returns the actual solution instead only its value.You may also want to get a copy of specific chromosome for a given population using get_chromosome.get_chromosome(brkga_data::BrkgaData,\n               population_index::Int64,\n               position::Int64)::Array{Float64, 1}For example, you can get the 3rd best chromosome from the 2nd population usingthird_best = get_chromosome(brkga_data, 2, 3)Now, suppose you get such chromosome or chromosomes and apply a quick local search procedure on them. It may be useful to reinsert such new solutions in the BRKGA population for the next evolutionary cycles. You can do that using inject_chromosome!.inject_chromosome!(brkga_data::BrkgaData,\n                   chromosome::Array{Float64, 1},\n                   population_index::Int64,\n                   position::Int64,\n                   fitness::Float64 = Inf)Note that the chromosome is put in a specific position of a given population. If you do not provide the fitness, inject_chromosome! will decode the injected chromosome. For example, the following code injects a random chromosome keys into the population #1 in the last position (population_size), i.e., it will replace the worst solution:keys = sort(rand(instance.num_nodes))\ninject_chromosome!(brkga_data, keys, 1, brkga_data.params.population_size)"
},

{
    "location": "guide/#Implicit-Path-Relink-1",
    "page": "Guide / Tutorial",
    "title": "Implicit Path Relink",
    "category": "section",
    "text": "The Implicit Path Relinking (IPR) is a nice addition to the standard BRKGA framework, and it provides an excellent way to create hybrid heuristics and push the optimization further. The good thing about IPR is that you do not need to worry about the path relink implementation, which can be long and tedious if done by hand or customized per problem.BrkgaMpIpr provides a friendly interface to use IPR directly from the BRKGA population, and you only must provide a few functions and arguments to have a Path Relink algorithm ready to go. This is the main signature of path_relink!path_relink!(brkga_data::BrkgaData,\n             pr_type::PathRelinkingType,\n             pr_selection::PathRelinkingSelection,\n             compute_distance::Function,\n             affect_solution::Function,\n             number_pairs::Int64,\n             minimum_distance::Float64,\n             block_size::Int64 = 1,\n             max_time::Float64 = 0.0,\n             percentage::Float64 = 1.0\n)::PathRelinkingResultThe first argument is the BrkgaData as usual. The 2nd argument defines the type of implicit path relink to be performed (PathRelinkingType). The DIRECT path relink exchanges the keys of two chromosomes directly, and it is usually more suitable to or threshold representations, i.e., where the key values are used to some kind of discretization, such as \" if x < 0.5, then 0, otherwise 1.\" The PERMUTATION path relink switches the order of a key according to its position in the other chromosome. Usually, this kind of path relink is more suitable to permutation representations, where the chromosome induces an order or permutation. For example, chromosome [0.4, 0.7, 0.1] may induce the increasing order (3, 1, 2). More details about threshold and permutation representations in this paper.PathRelinkingSelection defines how the algorithm picks the chromosomes for relinking. BESTSOLUTION selects, in the order, the best solution of each population. RANDOMELITE chooses uniformly random solutions from the elite sets.The next argument is a function to compute the distance between two chromosomes such signature must becompute_distance(vector1::Array{Float64, 1},\n                 vector2::Array{Float64, 1})::Float64If the value returned by compute_distance() is greater than or equal to minimum_distance, the algorithm will perform the path relink between the two chromosomes. Otherwise, it will look for another pair of chromosomes. The algorithm will try number_pairs chromosomes before gives up. In the presence of multiple populations, the path relinking is performed between elite chromosomes from different populations, in a circular fashion. For example, suppose we have 3 populations. The framework performs 3 path relinkings: the first between individuals from populations 1 and 2, the second between populations 2 and 3, and the third between populations 3 and 1. In the case of just one population, both base and guiding individuals are sampled from the elite set of that population.Note that in traditional path relink algorithms, compute_distance() depends on the problem structure. On IPR, you can use a generic distance function, or provide one that incorporates more knowledge about the problem. BrkgaMpIpr provides a function to compute the (modified) Hamming distance for threshold representations (hamming_distance), and a function that computes the Kendall Tau distance distance for permutation representations (kendall_tau_distance). Again, details about threshold and permutation representations in this paper.As a simple example, suppose you are using a threshold representation where each chromosome key can represent one of 3 different values (a ternary threshold representation). So, one possible way to compute the distance between two chromosomes can be:function value(key::Float64)::Float64\n    return key < 0.33 ? 0.0 : (key < 0.66 ? 1.0 : 2.0)\nend\n\nfunction compute_distance(vector1::Array{Float64, 1},\n                          vector2::Array{Float64, 1})::Float64\n    total = 0.0\n    for i in 1:length(vector1)\n        total += abs(value(vector1[i]) - value(vector2[i]))\n    end\n    return total\nendTo avoid changes that do not lead to new solutions, we must verify if such key exchanges affect the solution. For that, we must pass a function with the signature:affect_solution(block1::SubArray{Float64, 1},\n                block2::SubArray{Float64, 1})::Boolaffect_solution two gets partial chromosomes/block of genes block1 and block2 and checks whether changing the keys from block1 to block2 affects the solution. For instance, suppose that the alleles/keys are used as threshold such that values > 0.5 activate a feature. Suppose we have block1 = [0.3, 0.4, 0.1] and block2 = [0.4, 0.1, 0.2]. Since all values are below 0.5, changing the keys from block1 to block2 do not change the solution, and therefore, we can drop such change (and subsequently decoding). The blocks can hold only one key/allele, sequential key blocks, or even the whole chromosome. Note that affect_solution is crucial to the IPR performance since this function helps to avoid exploring regions already surveyed. Also, note that affect_solution can incorporate some problem knowledge.warning: Warning\nThe current implementation of permutation path relink does not make use of affect_solution. However, path_relink! requires the function. You can use the simple lambda function for this one:(x, y) -> trueblock_size defines the number of keys / size of the chromosome block to be exchanged during the direct path relink. This parameter is also critical for IPR performance since it avoids too many exchanges during the path building. Usually, we can compute this number based on the size of the chromosome by some factor (alpha_block_size in the configuration file), chosen by you. Again, details here.note: Note\nExperiments have shown that a good choice isblock_size = alpha_block_size times sqrtsizeofchromosomeThe last two parameters are stopping criteria. The algorithm stops either when max_time seconds is reached or percentage% of the path is built.warning: Warning\nIPR is a very time-intensive process. You must set the stopping criteria accordingly.Let\'s see the example on main_complete.jl. Remember, since we are solving the TSP, we want to use the permutation-based IPR, and the Kendall Tau distance functions.result = path_relink!(\n    brkga_data,\n    brkga_params.pr_type,\n    brkga_params.pr_selection,\n    kendall_tau_distance,\n    affect_solution_kendall_tau,\n    brkga_params.pr_number_pairs,\n    brkga_params.pr_minimum_distance,\n    1, #block_size doesn\'t not matter for permutation.\n    maximum_time - (time() - start_time),\n    brkga_params.pr_percentage\n)Note that most parameters come from BrkgaParams. The maximum IPR time is set to the remaining time for optimization (global maximum_time minus the elapsed time time() - start_time.path_relink! returns a PathRelinkingResult object which defines the status of the IPR optimization. These status are described on PathRelinkingResult.note: Note\nThe TOO_HOMOGENEOUS status is directly linked to the chosen distance function and minimum distance. If the minimum distance is too large, IPR may not be able to find a pair of chromosomes far enough for path relink.If the found solution is the best solution found so far, IPR replaces the worst solution by it. Otherwise, IPR computes the distance between the found solution and all other solutions in the elite set, and replaces the worst solution by it if and only if the found solution is, at least, minimum_distance from all them."
},

{
    "location": "guide/#Important-notes-about-IPR-1",
    "page": "Guide / Tutorial",
    "title": "Important notes about IPR",
    "category": "section",
    "text": "The API will call decode!() function always with writeback = false. The reason is that if the decoder rewrites the chromosome, the path between solutions is lost and inadvertent results may come up. Note that at the end of the path relinking, the method calls the decoder with writeback = true in the best chromosome found to guarantee that this chromosome is re-written to reflect the best solution found.warning: Warning\nMake sure your decoder does not rewrite the chromosome when called with the argument writeback = false.BrkgaMpIpr path_relink! implementation is multi-threaded. Instead of to build and decode each chromosome one at a time, the method builds a list of candidates, altering the alleles/keys according to the guide solution, and then decode all candidates in parallel. Note that O(chromosome_size^2  block_size) additional memory is necessary to build the candidates, which can be costly if the chromosome_size is very large.warning: Warning\nAs it is in evolve!(), the decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing)."
},

{
    "location": "guide/#Shaking-and-Resetting-1",
    "page": "Guide / Tutorial",
    "title": "Shaking and Resetting",
    "category": "section",
    "text": "Sometimes, BRKGA gets stuck, converging to local maxima/minima, for several iterations. When such a situation happens, it is a good idea to perturb the population, or even restart from a new one completely new. BrkgaMpIpr offers shake! function, an improved variation of the original version proposed in this paper.shake!(brkga_data::BrkgaData,\n       intensity::Int64,\n       shaking_type::ShakingType,\n       population_index::Int64 = Inf64)shake! function gets an intensity parameter that measures how many times the perturbation is applied on the elite set for a given population_index (if not given, all populations are shaken). This method offers two generic/implicit ShakingTypes. With CHANGE, direct modifications are done in the keys/alleles. This kind of shaking is recommended when the chromosome uses direct or threshold representations. SWAP exchanges keys/alleles inducing new permutations. For representational definitions, please read this paper. For instance, the following code shakes all populations using 10 swap moves.shake!(brkga_data, 10, SWAP)Sometimes, even shaking the populations does not help to escape from local maxima/minima. So, we need a drastic measure, restarting from scratch the role population. This can be easily accomplished with reset!.reset!(brkga_data)note: Note\nWhen using reset!, all warm-start solutions provided by set_initial_population! are discarded. You may use inject_chromosome! to insert those solutions again.warning: Warning\nAgain, the decoding of each chromosome is done in parallel if multi-thread is enabled. Therefore, we must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using a single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing)."
},

{
    "location": "guide/#Multi-population-and-migration-1",
    "page": "Guide / Tutorial",
    "title": "Multi-population and migration",
    "category": "section",
    "text": "Multi-population or island model was introduced in genetic algorithms in this paper. The idea is to evolve parallel and independent populations and, once a while, exchange individuals among these populations. In several scenarios, this approach is very beneficial for optimization.BrkgaMpIpr is implemented using such island idea from the core. If you read the guide until here, you may notice that several methods take into account multiple populations. To use multiple populations, you must set BrkgaParams.num_independent_populations with 2 ou more populations, and build BrkgaData from such parameters.The immigration process is implemented byexchange_elite!(brkga_data::BrkgaData, num_immigrants::Int64)exchange_elite! copies num_immigrants from one population to another, replacing the worst num_immigrants individuals from the recipient population. Note that the migration is done for all pairs of populations. For instance, the following code exchanges 3 best individuals from each population:exchange_elite!(brkga_data, 3)"
},

{
    "location": "guide/#Simulating-the-standard-BRKGA-1",
    "page": "Guide / Tutorial",
    "title": "Simulating the standard BRKGA",
    "category": "section",
    "text": "Sometimes, it is a good idea to test how the standard BRKGA algorithm performs for a problem. You can use BrkgaMpIpr framework to quickly implement and test a standard BRKGA.First, you must guarantee that, during the crossover, the algorithm chooses only one elite individual and only one non-elite individual. This is easily accomplished setting num_elite_parents = 1 and total_parents = 2. Then, you must set up a bias function that ranks the elite and no-elite individual according to the original BRKGA bias parameter rho (rho).You can use set_bias_custom_function! for that task. The given function receives the index of the chromosome and returns a ranking for it. Such ranking is used in the roulette method to choose the individual from which each allele comes to build the new chromosome. Since we have one two individuals for crossover in the standard BRKGA, the bias function must return the probability to one or other individual. In the following code, we do that with a simple if...else lambda function.# create brkga_params by hand or reading from a file,\n# then set the following by hand.\nbrkga_params.num_elite_parents = 1\nbrkga_params.total_parents = 2\n\nrho = 0.75\nset_bias_custom_function!(brkga_data, x -> x == 1 ? rho : 1.0 - rho)\ninitialize!(brkga_data)Here, we first set the num_elite_parents = 1 and total_parents = 2 as explained before. Following, we set a variable rho = 0.75. This is the rho from standard BRKGA, and you may set it as you wish. Then, we set the bias function as a very simple lambda function:x -> x == 1 ? rho : 1.0 - rhoSo, if the index of the chromosome is 1 (elite individual), it gets a 0.75 rank/probability. If the index is 2 (non-elite individual), the chromosome gets 0.25 rank/probability.note: Note\nAll these operations must be done before calling initialize!."
},

{
    "location": "guide/#Reading-and-writing-parameters-1",
    "page": "Guide / Tutorial",
    "title": "Reading and writing parameters",
    "category": "section",
    "text": "Although we can build a BrkgaData by set up a BrkgaParams object manually, the easiest way to do so is to read such parameters from a configuration file. For this, we can use load_configuration that reads a simple plain text file and returns a tuple of BrkgaParams and ExternalControlParams objects. For instance,brkga_params, control_params = load_configuration(\"tuned_conf.txt\")The configuration file must be plain text such that contains pairs of parameter name and value. This file must list all fields from BrkgaParams and ExternalControlParams, even though you do not use each one. In examples folder we have config.conf that looks like this:population_size 2000\nelite_percentage 0.30\nmutants_percentage 0.15\nnum_elite_parents 2\ntotal_parents 3\nbias_type LOGINVERSE\nnum_independent_populations 3\npr_number_pairs 0\npr_minimum_distance 0.15\npr_type PERMUTATION\npr_selection BESTSOLUTION\nalpha_block_size 1.0\npr_percentage 1.0\nexchange_interval 200\nnum_exchange_indivuduals 2\nreset_interval 600It does not matter whether we use lower or upper cases. Blank lines and lines starting with # are ignored. The order of the parameters should not matter either. And, finally, this file should be readble for both C++ and Julia framework versions.In some cases, you define some of the parameters at the running time, and you may want to save them for debug or posterior use. To do so, you can use write_configuration, call upon a BrkgaParams object or BrkgaData object. For example,write_configuration(\"my_new_parameters.conf\", brkga_params, external_params)\n# or\nwrite_configuration(\"crazy_parameters.txt\", brkga_data)note: Note\nwrite_configuration rewrites the given file. So, watch out to not lose previous configurations."
},

{
    "location": "guide/#(Probable-Valuable)-Tips-1",
    "page": "Guide / Tutorial",
    "title": "(Probable Valuable) Tips",
    "category": "section",
    "text": ""
},

{
    "location": "guide/#Algorithm-warmup-1",
    "page": "Guide / Tutorial",
    "title": "Algorithm warmup",
    "category": "section",
    "text": "When using Julia code, it is an excellent idea to dry-run all functions you may use and, mainly, want to time. The reason is that Julia uses lazy evaluation when live-compiling the code, i.e., it compiles as it goes. Another advantage is the memory location effects of our data (principle of locality), that can be brought closer to the processor (L2/L3 caches) during the running. Obliviously, this depends on how you implement and use your data structures.In main_complete.jl, we have the following piece of code to warmup mainly the decoder and other functions. Note that we just deep-copy brkga_data, and then, we may lose the principle of locality.bogus_data = deepcopy(brkga_data)\nevolve!(bogus_data, 2)\npath_relink!(bogus_data, brkga_params.pr_type, brkga_params.pr_selection,\n             (x, y) -> 1.0, (x, y) -> true, 0, 0.5, 1, 10.0, 1.0)\nget_best_fitness(brkga_data)\nget_best_chromosome(brkga_data)\nbogus_data = nothing"
},

{
    "location": "guide/#Complex-decoders-and-timing-1",
    "page": "Guide / Tutorial",
    "title": "Complex decoders and timing",
    "category": "section",
    "text": "Some problems require complex decoders while for others, the decoder contains local search procedures, that can be time-consuming. In general, the decoding is the most time-expensive component of a BRKGA algorithm, and it may skew some stopping criteria based on running time. Therefore, if your decoder is time-consuming, it is a good idea to implement a timer or chronometer kind of thing inside the decoder.Testing for stopping time uses several CPU cycles, and you need to be careful when/where to test it, otherwise, you spend all the optimization time doing system calls to the clock.IMHO, the most effective way to do it is to test time at the very end of the decoding. If the current time is larger than the maximum time allowed, simple return Inf or -Inf according to your optimization direction. In this way, we make the solution invalid since it violates the maximum time allowed. The BRKGA framework takes care of the rest."
},

{
    "location": "guide/#Multi-threading-1",
    "page": "Guide / Tutorial",
    "title": "Multi-threading",
    "category": "section",
    "text": "Since Moore\'s law is not holding its status anymore, we, simple mortals, must appeal to the wonders of multi-threading. This paradigm can be tricky to code, and Amdahl\'s law plays against us. Several genetic algorithms, and in particular, BRKGA, can use parallel solution evaluation (or decoding), which makes the use of multi-threading relatively straightforward. BrkgaMpIpr is not different, and it uses Julia multi-threading capabilities to do so.First, as commented several times in this guide, the decoder must be THREAD-SAFE. So, each thread must have its own read/write data structures and may share other read-only data. The simplest way to do it is to create those structures inside the decoder (like most people do). But be aware, this strategy slows down the algorithm significantly depending on the size and format of the structures, and I do not recommend it.IMHO, the best way to do that is to preallocate the data structure per thread (using Threads.nthreads()), and pass them to the decoder through the problem instance. Then, inside the decoder, you can use Threads.threadid() and recover the memory you want to use.Let\'s see a simple example considering the TSP example. tsp_decode! uses a single array to create the permutation of nodes. Let\'s pre-allocate its memory per thread. So, in TSP_Instance, we pre-allocate copies of such array, one for each thread:using Base.Threads\n\n# Declare as a type to make the code shorter and more readable.\nPermutationArray = Array{Tuple{Float64, Int64}, 1}\n\nstruct TSP_Instance <: AbstractInstance\n    num_nodes::Int64\n    distances::Array{Float64}\n\n    # Permutations arrays per thread, to be pre-allocated.\n    permutation_per_thread::Array{PermutationArray, 1}\n\n    function TSP_Instance(filename::String)\n        #... Code for loading here\n\n        # Allocate the main array to create references for other arrays per thread.\n        permutation_per_thread = Array{PermutationArray, 1}(undef, nthreads())\n\n        # Pre-allocate the permutation arrays, one for each thread.\n        for i in 1:nthreads()\n            permutation_per_thread[i] = PermutationArray(undef, num_nodes)\n        end\n\n        new(num_nodes, distances, permutation_per_thread)\n    end\nendThen, in tsp_decode!, we simply refer to the array according to the local thread ID:function tsp_decode!(chromosome::Array{Float64}, instance::TSP_Instance,\n                     rewrite::Bool = true)::Float64\n    permutation = instance.permutation_per_thread[threadid()]\n    #...\n    #...\nendNote that to pre-allocate decoding structures inside the object holding the instance is not the most elegant and decoupled code we can write. However, to decouple the decoding data from the instance data requires that we pass another data object to the decoder. To do this explicitly, we may get an embroidered API. We could do it implicitly, by creating a Singleton object to hold the decoding data. However, this also reduces (a lot) the clarity and objectivity of the code. In C++ code, this is much easier accomplished by creating a Decoder object that can hold data members as much as methods. Therefore, when creating the Decoder object, we can pre-allocate all data structures we need.note: Note\nPre-allocation and multi-threading only make sense for large data structures and time-consuming decoders. Otherwise, the code spends too much time on context switching and system calls."
},

{
    "location": "license/#",
    "page": "License",
    "title": "License",
    "category": "page",
    "text": ""
},

{
    "location": "license/#BRKGA-MP-IPR-License-1",
    "page": "License",
    "title": "BRKGA-MP-IPR License",
    "category": "section",
    "text": "Copyright (c) 2019, Carlos Eduardo de Andrade. All other rights reserved.Redistribution and use in source and binary forms, with or without modification, are only permitted provided that the following conditions are met:Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.\nRedistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.\nNeither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.\nAll publications, softwares, or any other materials mentioning features or use of this software and/or the data used to test it must cite explicitly the following article:C.E. Andrade. R.F. Toso, J.F. Gonçalves, M.G.C. Resende. The Multi-Parent Biased Random-key Genetic Algorithm with Implicit Path Relinking and its real-world applications.  European Journal of Operational Research, volume XX, issue X, pages xx-xx, 2019. DOI to be determinedTHIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
},

{
    "location": "api/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "api/#API-documentation-1",
    "page": "Library",
    "title": "API documentation",
    "category": "section",
    "text": ""
},

{
    "location": "api/#BrkgaMpIpr.BiasFunction",
    "page": "Library",
    "title": "BrkgaMpIpr.BiasFunction",
    "category": "type",
    "text": "@enum BiasFunction\n\nSpecifies a bias function when choosing parents to mating. This function substitutes the rho (rho) parameter from the original BRKGA. For a given rank r, we have the following functions:\n\nCONSTANT: 1 / number of parents for mating (all individuals have the             same probability)\nCUBIC: r^-3\nEXPONENTIAL: ϵ^-r\nLINEAR: 1  r\nLOGINVERSE: 1  log(r + 1)\nQUADRATIC: r^-2\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.PathRelinkingResult",
    "page": "Library",
    "title": "BrkgaMpIpr.PathRelinkingResult",
    "category": "type",
    "text": "@enum PathRelinkingResult\n\nSpecifies the result type/status of path relink procedure:\n\nTOO_HOMOGENEOUS: the chromosomes among the populations are too homogeneous                    and the path relink will not generate improveded                    solutions.\nNO_IMPROVEMENT: path relink was done but no improveded solution was found.\nELITE_IMPROVEMENT: an improved solution among the elite set was found,                      but the best solution was not improved.\nBEST_IMPROVEMENT: the best solution was improved.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.PathRelinkingSelection",
    "page": "Library",
    "title": "BrkgaMpIpr.PathRelinkingSelection",
    "category": "type",
    "text": "@enum PathRelinkingSelection\n\nSpecifies which individuals used to build the path:\n\nBESTSOLUTION: selects, in the order, the best solution of each population.\nRANDOMELITE: chooses uniformly random solutions from the elite sets.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.PathRelinkingType",
    "page": "Library",
    "title": "BrkgaMpIpr.PathRelinkingType",
    "category": "type",
    "text": "@enum PathRelinkingType\n\nSpecifies type of path relinking:\n\nDIRECT: changes each key for the correspondent in the other chromosome.\nPERMUTATION: switches the order of a key for that in the other chromosome.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.Sense",
    "page": "Library",
    "title": "BrkgaMpIpr.Sense",
    "category": "type",
    "text": "@enum Sense\n\nTells the algorithm either to MINIMIZE or MAXIMIZE the objective function.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.ShakingType",
    "page": "Library",
    "title": "BrkgaMpIpr.ShakingType",
    "category": "type",
    "text": "@enum ShakingType\n\nSpecifies the type of shaking to be performed.\n\nCHANGE: applies the following perturbations:\nInverts the value of a random chosen, i.e., from value to 1 - value;\nAssigns a random value to a random key.\nSWAP: applies two swap perturbations:\nSwaps the values of a randomly chosen key i and its neighbor i + 1;\nSwaps values of two randomly chosen keys.\n\n\n\n\n\n"
},

{
    "location": "api/#Enumerations-1",
    "page": "Library",
    "title": "Enumerations",
    "category": "section",
    "text": "BiasFunction\nPathRelinkingResult\nPathRelinkingSelection\nPathRelinkingType\nSense\nShakingType"
},

{
    "location": "api/#BrkgaMpIpr.AbstractInstance",
    "page": "Library",
    "title": "BrkgaMpIpr.AbstractInstance",
    "category": "type",
    "text": "abstract type AbstractInstance\n\nThe required interface for external data to be provided to the decoder. The problem definitions and data must be a subtype of AbstractInstance. For example,\n\nmutable struct TSPInstance <: AbstractInstance\n    num_cities::Int64\n    distances::Array{Float64}\nend\n\nrepresents an instance type for the Traveling Salesman problem which defines the number fo cities and a matrix of distances between them.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.BrkgaData",
    "page": "Library",
    "title": "BrkgaMpIpr.BrkgaData",
    "category": "type",
    "text": "mutable struct BrkgaData\n\nRepresents the internal state of the BRKGA-MP-IPR algorithm.\n\nThis structure has no direct constructor and must be built using build_brkga functions. You can create multiple BrkgaData representing different states of the algorithm, and use them independently.\n\nwarning: Warning\nThis structure is NOT INTENDED to be used outside BRKGA functions. Ad hoc changes may lead to inadvertent results.\n\nFields\n\nopt_sense\nOptimization sense (minimization = 0, maximization = 1).\n\nchromosome_size\nNumber of genes in the chromosome [> 0].\n\nparams\nBRKGA parameters for evolution.\n\nelite_size\nNumber of elite items in the population [> 0].\n\nnum_mutants\nNumber of mutants introduced at each generation into the population [> 0].\n\nevolutionary_mechanism_on\nIf false, no evolution is performed but only chromosome decoding. Very useful to emulate a multi-start algorithm.\n\nproblem_instance\n(Internal data) The problem instance used by the decode! function to map a chromosome to a problem solution. Since decode! should not change this data, this attribute can be considered constant.\n\ndecode!\n(Internal data) This is the main decode function called during the evolutionary process and in the path relink. It must have the following signature:\ndecode!(chromosome::Array{Float64, 1},\n        problem_instance::AbstractInstance,\n        rewrite::Bool = true)::Float64\nNote that if rewrite == false, decode! must not change chromosome. IPR routines requires decode! to not change chromosome.\n\nrng\n(Internal data) The internal random generator number. DON\'T USE FOR ANYTHING OUTSIDE. If you need a RNG, create a new generator.\n\nprevious\n(Internal data) Previous population.\n\ncurrent\n(Internal data) Current population.\n\nbias_function\n(Internal data) A unary non-increasing function such that bias_function(::Int64 > 0)::Float64\n\ntotal_bias_weight\n(Internal data) Holds the sum of the results of each raking given a bias function. This value is needed to normalization.\n\nshuffled_individuals\n(Internal data) Used to shuffled individual/chromosome indices during the mate.\n\nparents_ordered\n(Internal data) Defines the order of parents during the mating.\n\ninitialized\n(Internal data) Indicates if the algorithm was proper initialized.\n\nreset_phase\n(Internal data) Indicates if the algorithm have been reset.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.BrkgaParams",
    "page": "Library",
    "title": "BrkgaMpIpr.BrkgaParams",
    "category": "type",
    "text": "mutable struct BrkgaParams\n\nRepresents the BRKGA and IPR hyper-parameters. You can load these parameters from a configuration file using load_configuration and build_brkga, and write them using write_configuration.\n\nFields\n\npopulation_size\nNumber of elements in the population [> 0].\n\nelite_percentage\nPercentage of individuals to become the elite set (0, 1].\n\nmutants_percentage\nPercentage of mutants to be inserted in the population\n\nnum_elite_parents\nNumber of elite parents for mating [> 0].\n\ntotal_parents\nNumber of total parents for mating [> 0].\n\nbias_type\nType of bias that will be used.\n\nnum_independent_populations\nNumber of independent parallel populations.\n\npr_number_pairs\nNumber of pairs of chromosomes to be tested to path relinking.\n\npr_minimum_distance\nMininum distance between chromosomes selected to path-relinking.\n\npr_type\nPath relinking type.\n\npr_selection\nIndividual selection to path-relinking.\n\nalpha_block_size\nDefines the block size based on the size of the population.\n\npr_percentage\nPercentage / path size to be computed. Value in (0, 1].\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.ExternalControlParams",
    "page": "Library",
    "title": "BrkgaMpIpr.ExternalControlParams",
    "category": "type",
    "text": "mutable struct ExternalControlParams\n\nRepresents additional control parameters that can be used outside this framework. You can load these parameters from a configuration file using load_configuration and build_brkga, and write them using write_configuration.\n\nFields\n\nexchange_interval\nInterval at which elite chromosomes are exchanged (0 means no exchange) [> 0].\n\nnum_exchange_indivuduals\nNumber of elite chromosomes exchanged from each population [> 0].\n\nreset_interval\nInterval at which the populations are reset (0 means no reset) [> 0].\n\n\n\n\n\n"
},

{
    "location": "api/#Types-1",
    "page": "Library",
    "title": "Types",
    "category": "section",
    "text": "AbstractInstance\nBrkgaData\nBrkgaParams\nExternalControlParams"
},

{
    "location": "api/#Base.parse",
    "page": "Library",
    "title": "Base.parse",
    "category": "function",
    "text": "parse(::Type{BiasFunction}, value::String)::BiasFunction\n\nParse value returning a valid BiasFunction enumeration.\n\nThrows\n\nArgumentError: in case the bias description does not match.\n\n\n\n\n\nparse(::Type{PathRelinkingType}, value::String)::PathRelinkingType\n\nParse value returning a valid PathRelinkingType enumeration.\n\nThrows\n\nArgumentError: in case the type description does not match.\n\n\n\n\n\nparse(::Type{PathRelinkingSelection}, value::String)::PathRelinkingSelection\n\nParse value returning a valid PathRelinkingSelection enumeration.\n\nThrows\n\nArgumentError: in case the selection description does not match.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.load_configuration",
    "page": "Library",
    "title": "BrkgaMpIpr.load_configuration",
    "category": "function",
    "text": "load_configuration(configuration_file::String)::\n        Tuple{BrkgaParams, ExternalControlParams}\n\nLoad the parameters from filename returning them as a tuple.\n\nThrows\n\nLoadError: in cases of the file is an invalid configuration file, parameters are missing, or parameters are ill-formatted.\nSystemError: in case the configuration files cannot be openned.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.write_configuration",
    "page": "Library",
    "title": "BrkgaMpIpr.write_configuration",
    "category": "function",
    "text": "function write_configuration(filename::String, brkga_params::BrkgaParams,\n                             external_params::ExternalControlParams)\n\nWrite brkga_params and external_params into filename.\n\nThrows\n\nSystemError: in case the configuration files cannot be openned.\n\n\n\n\n\nfunction write_configuration(filename::String, brkga_data::BrkgaData,\n                             external_params::ExternalControlParams =\n                                                    ExternalControlParams())\n\nWrite the parameters from brkga_data.params and external_params into filename.\n\nThrows\n\nSystemError: in case the configuration files cannot be openned.\n\n\n\n\n\n"
},

{
    "location": "api/#I/O-functions-1",
    "page": "Library",
    "title": "I/O functions",
    "category": "section",
    "text": "parse\nload_configuration\nwrite_configuration"
},

{
    "location": "api/#BrkgaMpIpr.build_brkga",
    "page": "Library",
    "title": "BrkgaMpIpr.build_brkga",
    "category": "function",
    "text": "build_brkga(problem_instance::AbstractInstance, decode_function!::Function,\n            opt_sense::Sense, seed::Int64, chromosome_size::Int64,\n            brkga_params::BrkgaParams, evolutionary_mechanism_on::Bool = true,\n)::BrkgaData\n\nBuild a BrkgaData object to be used in the evolutionary and path relink procedures. Such data structure should not be changed outside the BrkgaMpIpr functions. This version accepts all control arguments, and it is handy for tuning purposes.\n\nArguments\n\nproblem_instance::AbstractInstance: an instance to the problem to be solved.\ndecode_function!::Function: the decode funtion used to map chromosomes to solutions. It must have the following signature:\ndecode!(chromosome::Array{Float64, 1},\n        problem_instance::AbstractInstance,\n        rewrite::Bool = true)::Float64\nNote that if rewrite == false, decode! cannot modify chromosome.\nopt_sense::Sense: the optimization sense ( maximization or minimization).\nseed::Int64: seed for the random number generator.\nchromosome_size::Int64: number of genes in each chromosome.\nbrkga_params::BrkgaParams: BRKGA and IPR parameters object loaded from a configuration file or manually created. All the data is deep-copied.\nevolutionary_mechanism_on::Bool = true: if false, no evolution is performed but only chromosome decoding. On each generation, all population is replaced excluding the best chromosome. Very useful to emulate a multi-start algorithm.\n\nThrows\n\nArgumentError: in several cases where the arguments or a combination of them are invalid.\n\n\n\n\n\nbuild_brkga(problem_instance, decode_function!, opt_sense, seed,\n    chromosome_size, configuration_file,\n    evolutionary_mechanism_on)::Tuple{BrkgaData, ExternalControlParams}\n\nBuild a BrkgaData object to be used in the evolutionary and path relink procedures, and a ExternalControlParams that holds additional control parameters. Note that BrkgaData should not be changed outside the BrkgaMpIpr functions. This version reads most of the parameters from a configuration file.\n\nArguments\n\nproblem_instance::AbstractInstance: an instance to the problem to be solved.\ndecode_function!::Function: the decode funtion used to map chromosomes to solutions. It must have the following signature:\ndecode!(chromosome::Array{Float64, 1},\n        problem_instance::AbstractInstance,\n        rewrite::Bool = true)::Float64\nNote that if rewrite == false, decode! cannot modify chromosome.\nopt_sense::Sense: the optimization sense ( maximization or minimization).\nseed::Int64: seed for the random number generator.\nchromosome_size::Int64: number of genes in each chromosome..\nconfiguration_file::String:  text file with the BRKGA parameters. An example can be found at <a href=\"example.conf\">example.conf</a>. Note that the content after \"#\" is considered comments and it is ignored.\nevolutionary_mechanism_on::Bool = true: if false, no evolution is performed but only chromosome decoding. On each generation, all population is replaced excluding the best chromosome. Very useful to emulate a multi-start algorithm.\n\nThrows\n\nLoadError: in cases of the file is an invalid configuration file, parameters are missing, or parameters are ill-formatted.\nSystemError: in case the configuration files cannot be openned.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.initialize!",
    "page": "Library",
    "title": "BrkgaMpIpr.initialize!",
    "category": "function",
    "text": "initialize!(brkga_data::BrkgaData)\n\nInitialize the populations and others data structures of the BRKGA. If an initial population is supplied, this method completes the remaining individuals, if they do not exist.\n\nwarning: Warning\nTHIS METHOD MUST BE CALLED BEFORE ANY OPTIMIZATION METHODS.\n\nThis method also performs the initial decoding of the chromosomes. Therefore, depending on the decoder implementation, this can take a while, and the user may want to time such procedure in his/her experiments.\n\nnote: Note\nAs it is in evolve!, the decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using a single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing).\n\nThrows\n\nErrorException: if bias_function is not defined previously.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.set_bias_custom_function!",
    "page": "Library",
    "title": "BrkgaMpIpr.set_bias_custom_function!",
    "category": "function",
    "text": "set_bias_custom_function!(brkga_data::BrkgaData, bias_function::Function)\n\nSet a new bias function to be used to rank the chromosomes during the mating. It must be a positive non-decreasing function returning a Float64, i.e., f(::Int64)::Float64 such that f(i) ≥ 0 and f(i) ≥ f(i+1) for i ∈ [1..total_parents]. For instance, the following sets an inverse quadratic function:\n\nset_bias_custom_function!(brkga_data, x -> 1.0 / (x * x))\n\nThrows\n\nArgumentError: in case the function is not a non-decreasing positive function.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.set_initial_population!",
    "page": "Library",
    "title": "BrkgaMpIpr.set_initial_population!",
    "category": "function",
    "text": "set_initial_population!(brkga_data::BrkgaData,\n                        chromosomes::Array{Array{Float64, 1}, 1})\n\nSet initial individuals into the poulation to work as warm-starters. Such individuals can be obtained from solutions of external procedures such as fast heuristics, other methaheuristics, or even relaxations from a mixed integer programming model that models the problem.\n\nAll given solutions are assigned to one population only. Therefore, the maximum number of solutions is the size of the populations.\n\nThrows\n\nArgumentError: if the number of given chromosomes is larger than the population size; if the sizes of the given chromosomes do not match with the required chromosome size.\n\n\n\n\n\n"
},

{
    "location": "api/#building_funcs-1",
    "page": "Library",
    "title": "Building functions",
    "category": "section",
    "text": "build_brkga\ninitialize!\nset_bias_custom_function!\nset_initial_population!"
},

{
    "location": "api/#BrkgaMpIpr.exchange_elite!",
    "page": "Library",
    "title": "BrkgaMpIpr.exchange_elite!",
    "category": "function",
    "text": "exchange_elite!(brkga_data::BrkgaData, num_immigrants::Int64)\n\nExchange elite-solutions between the populations. Given a population, the num_immigrants best solutions are copied to the neighbor populations, replacing their worth solutions. If there is only one population, nothing is done.\n\nThrows\n\nErrorException: if initialize! has not been called before.\nArgumentError: when num_immigrants < 1.\nArgumentError: num_immigrants ≥ ⌈population_size/num_independent_populations⌉ - 1.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.inject_chromosome!",
    "page": "Library",
    "title": "BrkgaMpIpr.inject_chromosome!",
    "category": "function",
    "text": "function inject_chromosome!(brkga_data::BrkgaData,\n                            chromosome::Array{Float64, 1},\n                            population_index::Int64,\n                            position::Int64,\n                            fitness::Float64 = Inf)\n\nInject chromosome and its fitness into population population_index in the position place. If fitness is not provided (fitness = Inf), the decoding is performed over chromosome. Once the chromosome is injected, the population is re-sorted according to the chromosomes\' fitness.\n\nThrows\n\nErrorException: if initialize! has not been called before.\nArgumentError: when population_index < 1 or population_index > num_independent_populations.\nArgumentError: when position < 1 or position > population_size.\nArgumentError: when lenght(chromosome) != chromosome_size.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.reset!",
    "page": "Library",
    "title": "BrkgaMpIpr.reset!",
    "category": "function",
    "text": "reset!(brkga_data::BrkgaData)\n\nReset all populations with brand new keys. All warm-start solutions provided by set_initial_population! are discarded. You may use inject_chromosome! to insert those solutions again.\n\nwarning: Warning\nAs it is in evolve!(), the decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing)\n\nThrows\n\nErrorException: if initialize! has not been called before.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.shake!",
    "page": "Library",
    "title": "BrkgaMpIpr.shake!",
    "category": "function",
    "text": "function shake!(brkga_data::BrkgaData, intensity::Int64,\n                shaking_type::ShakingType, population_index::Int64 = Inf64)\n\nPerform a shaking in the chosen population. The procedure applies changes (shaking) on elite chromosomes and fully reset the remaining population.\n\nArguments\n\nintensity::Int64: the intensity of the shaking (> 0);\nshaking_type::ShakingType: either CHANGE or SWAP moves;\npopulation_index::Int64: the index of the population to be shaken. If population_index > num_independent_populations, all populations are shaken.\n\nThrows\n\nErrorException: if initialize! has not been called before.\nArgumentError: when population_index < 1 or intensity < 1.\n\n\n\n\n\n"
},

{
    "location": "api/#Population-manipulation-functions-1",
    "page": "Library",
    "title": "Population manipulation functions",
    "category": "section",
    "text": "exchange_elite!\ninject_chromosome!\nreset!\nshake!"
},

{
    "location": "api/#BrkgaMpIpr.get_best_chromosome",
    "page": "Library",
    "title": "BrkgaMpIpr.get_best_chromosome",
    "category": "function",
    "text": "get_best_chromosome(brkga_data::BrkgaData)::Array{Float64, 1}\n\nReturn a copy of the best individual found so far among all populations.\n\nThrows\n\nErrorException: if initialize! has not been called before.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.get_best_fitness",
    "page": "Library",
    "title": "BrkgaMpIpr.get_best_fitness",
    "category": "function",
    "text": "get_best_fitness!(brkga_data::BrkgaData)::Float64\n\nReturn the fitness/value of the best individual found so far among all populations.\n\nThrows\n\nErrorException: if initialize! has not been called before.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.get_chromosome",
    "page": "Library",
    "title": "BrkgaMpIpr.get_chromosome",
    "category": "function",
    "text": "get_chromosome(brkga_data::BrkgaData, population_index::Int64,\n               position::Int64)::Array{Float64, 1}\n\nReturn a copy of the chromosome position in the population population_index.\n\nThrows\n\nErrorException: if initialize! has not been called before.\nArgumentError: when population_index < 1 or population_index > num_independent_populations.\nArgumentError: when position < 1 or position > population_size.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.get_current_population",
    "page": "Library",
    "title": "BrkgaMpIpr.get_current_population",
    "category": "function",
    "text": "get_current_population(brkga_data::BrkgaData,\n                       population_index::Int64)::Population\n\nReturn a reference for population population_index.\n\nnote: Note\nThis function is implemented for complaince with the C++ API. The user can access the population directly using brkga_data.current[population_index].\n\nwarning: Warning\nIT IS NOT ADIVISED TO CHANGE THE POPULATION DIRECTLY, since such changes can result in undefined behavior.\n\nThrows\n\nErrorException: if initialize! has not been called before.\nArgumentError: when population_index < 1 or population_index > num_independent_populations.\n\n\n\n\n\n"
},

{
    "location": "api/#Retrival-functions-1",
    "page": "Library",
    "title": "Retrival functions",
    "category": "section",
    "text": "get_best_chromosome\nget_best_fitness\nget_chromosome\nget_current_population"
},

{
    "location": "api/#BrkgaMpIpr.evolve!",
    "page": "Library",
    "title": "BrkgaMpIpr.evolve!",
    "category": "function",
    "text": "evolve!(brkga_data::BrkgaData, num_generations::Int64 = 1)\n\nEvolve the current populations following the guidelines of Multi-parent BRKGAs for num_generations generations.\n\nwarning: Warning\nThe decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing).\n\nThrows\n\nErrorException: if initialize!() was not called before.\nArgumentError: when num_generations < 1.\n\n\n\n\n\n"
},

{
    "location": "api/#Evolution-functions-1",
    "page": "Library",
    "title": "Evolution functions",
    "category": "section",
    "text": "evolve!"
},

{
    "location": "api/#BrkgaMpIpr.affect_solution_hamming_distance",
    "page": "Library",
    "title": "BrkgaMpIpr.affect_solution_hamming_distance",
    "category": "function",
    "text": "affect_solution_hamming_distance(block1::SubArray{Float64, 1},\n                                 block2::SubArray{Float64, 1},\n                                 threshold::Float64 = 0.5)::Bool\n\nReturn true the the changing of the blocks of keys block1 by the blocks of keys block2 affects the solution, based on the Hamming distance.\n\nnote: Note\nThis function may be more appropriated to threshold/direct chromosome representations.\n\nnote: Note\nblock1 and block2 must have the same size. No bounds checking is done due to performance reasons.\n\nnote: Note\nThis function is annotated with @inline due to performance reasons too.\n\nArguments\n\nblock1::SubArray{Float64, 1}: the first vector.\nblock2::SubArray{Float64, 1}: the second vector.\nthreshold::Float64 = 0.5: the threshold for binarization.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.affect_solution_kendall_tau",
    "page": "Library",
    "title": "BrkgaMpIpr.affect_solution_kendall_tau",
    "category": "function",
    "text": "affect_solution_kendall_tau(block1::SubArray{Float64, 1},\n                            block2::SubArray{Float64, 1})::Bool\n\nReturn true the the changing of the blocks of keys block1 by the blocks of keys block2 affects the solution, based on the Kendall Tau distance.\n\nnote: Note\nblock1 and block2 must have the same size. No bounds checking is done due to performance reasons.\n\nArguments\n\nblock1::SubArray{Float64, 1}: the first vector.\nblock2::SubArray{Float64, 1}: the second vector.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.hamming_distance",
    "page": "Library",
    "title": "BrkgaMpIpr.hamming_distance",
    "category": "function",
    "text": "hamming_distance(vector1::Array{Float64, 1}, vector2::Array{Float64, 1},\n                 threshold::Float64 = 0.5)::Float64\n\nCompute the Hamming distance between two vectors. It takes a threshold parameter to \"binarize\" the vectors. For instance, if threshold = 0.7, all values larger than or equal to 0.7 will be considerd 1.0, otherwise 0.0.\n\nnote: Note\nThis function may be more appropriated to threshold/direct chromosome representations.\n\nArguments\n\nvector1::Array{Float64, 1}: the first vector.\nvector2::Array{Float64, 1}: the second vector.\nthreshold::Float64 = 0.5: the threshold for binarization.\n\nThrows\n\nArgumentError: if vector1 and vector2 have different sizes.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.kendall_tau_distance",
    "page": "Library",
    "title": "BrkgaMpIpr.kendall_tau_distance",
    "category": "function",
    "text": "kendall_tau_distance(vector1::Array{Float64, 1},\n                     vector2::Array{Float64, 1};\n                     stop_immediately::Bool = false)::Float64\n\nkendall_tau_distance(vector1::SubArray{Float64, 1},\n                     vector2::SubArray{Float64, 1};\n                     stop_immediately::Bool = false)::Float64\n\nCompute the Kendall Tau distance between two vectors.\n\nnote: Note\nThis function may be more appropriated to permutation chromosome representations.\n\nArguments\n\nvector1::Array{Float64, 1}: the first vector.\nvector2::Array{Float64, 1}: the second vector.\nstop_immediately::Bool = false: if true, stop the computation immediately after find a difference.\n\nThrows\n\nArgumentError: if vector1 and vector2 have different sizes.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.path_relink!",
    "page": "Library",
    "title": "BrkgaMpIpr.path_relink!",
    "category": "function",
    "text": "function path_relink!(brkga_data::BrkgaData,\n    pr_type::PathRelinkingType,\n    pr_selection::PathRelinkingSelection,\n    compute_distance::Function,\n    affect_solution::Function,\n    number_pairs::Int64,\n    minimum_distance::Float64,\n    block_size::Int64,\n    max_time::Float64,\n    percentage::Float64\n)::PathRelinkingResult\n\nPerform path relinking between elite solutions that are, at least, a given minimum distance between themselves.\n\nIn the presence of multiple populations, the path relinking is performed between elite chromosomes from different populations, in a circular fashion. For example, suppose we have 3 populations. The framework performs 3 path relinkings: the first between individuals from populations 1 and 2, the second between populations 2 and 3, and the third between populations 3 and 1. In the case of just one population, both base and guiding individuals are sampled from the elite set of that population.\n\nNote that the algorithm tries to find a pair of base and guiding solutions with a minimum distance given by the distance function. If this is not possible, a new pair of solutions are sampled (without replacement) and tested against the distance. In case it is not possible to find such pairs for the given populations, the algorithm skips to the next pair of populations (in a circular fashion, as described above). Yet, if such pairs are not found in any case, the algorithm declares failure. This indicates that the populations are very homogeneous.\n\nIf the found solution is the best solution found so far, IPR replaces the worst solution by it. Otherwise, IPR computes the distance between the found solution and all other solutions in the elite set, and replaces the worst solution by it if and only if the found solution is, at least, minimum_distance from all them.\n\nThe API will call decode!() function always with writeback = false. The reason is that if the decoder rewrites the chromosome, the path between solutions is lost and inadvertent results may come up. Note that at the end of the path relinking, the method calls the decoder with writeback = true in the best chromosome found to guarantee that this chromosome is re-written to reflect the best solution found.\n\nThis method is a multi-thread implementation. Instead of to build and decode each chromosome one at a time, the method builds a list of candidates, altering the alleles/keys according to the guide solution, and then decode all candidates in parallel. Note that O(chromosome_size^2 / block_size) additional memory is necessary to build the candidates, which can be costly if the chromosome_size is very large.\n\nwarning: Warning\nAs it is in evolve!(), the decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing).\n\nArguments\n\nbrkga_data::BrkgaData: the BRKGA data.\npr_type::PathRelinkingType: type of path relinking to be performed. Either DIRECT or PERMUTATION-based.\npr_selection::PathRelinkingSelection: selection of which individuals use to path relinking. Either BESTSOLUTION or RANDOMELITE.\ncompute_distance::Function: the function used to compute the distance between two chromosomes. The function MUST HAVE the following signature\ncompute_distance(vector1::Array{Float64, 1},\n                 vector2::Array{Float64, 1})::Float64\naffect_solution::Function: function that takes two partial chromosomes / block of genes block1 and block2 and checks whether changing the keys from block1 to block2 affects the solution. For instance, suppose that the alleles/keys are used as threshold such that values > 0.5 activate a feature. Suppose we have block1 = [0.3, 0.4, 0.1] and block2 = [0.4, 0.1, 0.2]. Since all values are below 0.5, changing the keys from block1 to block2 do not change the solution, and therefore, we can drop such change (and subsequentely decoding). The blocks can hold only one key/allele, sequential key blocks, or even the whole chromosome. affect_solution takes two views/subarrays. The function MUST HAVE the following signature\naffect_solution(block1::SubArray{Float64, 1},\n                block2::SubArray{Float64, 1})::Bool\nnote: Note\nThis function depends on the problem structure and how the   keys/alleles are used\nnumber_pairs::Int64: number of chromosome pairs to be tested.  If number_pairs < 1, all pairs are tested.\nminimum_distance::Float64: minimum distance between two chromosomes computed by compute_distance.\nblock_size::Int64 = 1: number of alleles to be exchanged at once in each iteration. If one, the traditional path relinking is performed. It must be ≥ 1.\nmax_time::Float64 = 0: abort path-relinking when reach max_time. If max_time ≤ 0, no limit is imposed. Given in seconds.\npercentage::Float64 = 1.0: define the size, in percentage, of the path to build. Range [0, 1].\n\nReturns\n\nReturns PathRelinkingResult depending on the relink status.\n\nThrows\n\nErrorException: if initialize!() was not called before.\nArgumentError: when percentage < 1e-6 || percentage > 1.0 and block_size < 1.\n\n\n\n\n\n"
},

{
    "location": "api/#Path-relink-functions-1",
    "page": "Library",
    "title": "Path relink functions",
    "category": "section",
    "text": "affect_solution_hamming_distance\naffect_solution_kendall_tau\nhamming_distance\nkendall_tau_distance\npath_relink!"
},

{
    "location": "api/#Internals-1",
    "page": "Library",
    "title": "Internals",
    "category": "section",
    "text": "These types and functions are used internally in the framework. They are not meant to be used directly.CurrentModule = BrkgaMpIpr"
},

{
    "location": "api/#BrkgaMpIpr.DecodeStruct",
    "page": "Library",
    "title": "BrkgaMpIpr.DecodeStruct",
    "category": "type",
    "text": "mutable struct DecodeStruct\n\nHold the data structures used to build a candidate chromosome for parallel decoding on permutation-based path relink.\n\nwarning: Warning\nTHIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.Population",
    "page": "Library",
    "title": "BrkgaMpIpr.Population",
    "category": "type",
    "text": "mutable struct Population (internal BRKGA data struct)\n\nEncapsulates a population of chromosomes. Note that this struct is NOT meant to be used externally of this unit.\n\nFields\n\nchromosomes\nPopulation of chromosomes.\n\nfitness\nFitness of a each chromosome. Each pair represents the fitness and the chromosome index.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.Triple",
    "page": "Library",
    "title": "BrkgaMpIpr.Triple",
    "category": "type",
    "text": "mutable struct Triple\n\nHold the data structures used to build a candidate chromosome for parallel decoding on direct path relink.\n\nwarning: Warning\nTHIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.\n\n\n\n\n\n"
},

{
    "location": "api/#Types-2",
    "page": "Library",
    "title": "Types",
    "category": "section",
    "text": "DecodeStruct\nPopulation\nTriple"
},

{
    "location": "api/#Base.:|",
    "page": "Library",
    "title": "Base.:|",
    "category": "function",
    "text": "Base.:|(x::PathRelinkingResult,\n        y::PathRelinkingResult)::PathRelinkingResult\n\nPerform bitwise OR between two PathRelinkingResult returning the highest rank PathRelinkingResult.\n\nExamples\n\njulia> TOO_HOMOGENEOUS | NO_IMPROVEMENT\nNO_IMPROVEMENT::PathRelinkingResult = 1\n\njulia> NO_IMPROVEMENT | ELITE_IMPROVEMENT\nELITE_IMPROVEMENT::PathRelinkingResult = 3\n\njulia> ELITE_IMPROVEMENT | BEST_IMPROVEMENT\nBEST_IMPROVEMENT::PathRelinkingResult = 7\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.empty_function",
    "page": "Library",
    "title": "BrkgaMpIpr.empty_function",
    "category": "function",
    "text": "const empty_function() = nothing\n\nRepresent an empty function to be used as flag during data and bias function setups.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.find_block_range",
    "page": "Library",
    "title": "BrkgaMpIpr.find_block_range",
    "category": "function",
    "text": "find_block_range(block_number::Int64, block_size::Int64,\n                 max_end::Int64)::UnitRange{Int64}\n\nReturn a positive range for the given block_number with length block_size, limited to the max_end.\n\nnote: Note\nThis function only accept positive numbers, and all sanity check is disregarded due to performance reasons.\n\nwarning: Warning\nTHIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.swap!",
    "page": "Library",
    "title": "BrkgaMpIpr.swap!",
    "category": "function",
    "text": "swap!(x::Array{Any, 1}, pos1::Int64, pos2::Int64)\n\nSwap the value in position pos1 with the value in position pos2 in vector x.\n\nnote: Note\nThis function only accept positive numbers, and all sanity and bounds check is disregarded due to performance reasons.\n\nwarning: Warning\nTHIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.\n\n\n\n\n\n"
},

{
    "location": "api/#Minor-helper-functions-1",
    "page": "Library",
    "title": "Minor helper functions",
    "category": "section",
    "text": ":|\nempty_function\nfind_block_range\nswap!"
},

{
    "location": "api/#BrkgaMpIpr.evolve_population!",
    "page": "Library",
    "title": "BrkgaMpIpr.evolve_population!",
    "category": "function",
    "text": "evolve_population!(brkga_data::BrkgaData, population_index::Int64)\n\nEvolve the population population_index to the next.\n\nwarning: Warning\nThe decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing).\n\nThrows\n\nErrorException: if initialize!() was not called before.\nArgumentError: when population_index < 1 or population_index > num_independent_populations.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.direct_path_relink!",
    "page": "Library",
    "title": "BrkgaMpIpr.direct_path_relink!",
    "category": "function",
    "text": "function direct_path_relink!(brkga_data::BrkgaData,\n                             chromosome1::Array{Float64, 1},\n                             chromosome2::Array{Float64, 1},\n                             affect_solution::Function,\n                             block_size::Int64,\n                             max_time::Float64,\n                             percentage::Float64\n    )::Tuple{Float64, Array{Float64, 1}}\n\nPerform the direct path relinking, changing each allele or block of alleles of base chromosome for the correspondent one in the guide chromosome.\n\nThe API will call decode!() function always with writeback = false. The reason is that if the decoder rewrites the chromosome, the path between solutions is lost and inadvertent results may come up. Note that at the end of the path relinking, the method calls the decoder with writeback = true in the best chromosome found to guarantee that this chromosome is re-written to reflect the best solution found.\n\nThis method is a multi-thread implementation. Instead of to build and decode each chromosome one at a time, the method builds a list of candidates, altering the alleles/keys according to the guide solution, and then decode all candidates in parallel. Note that O(chromosome_size^2 / block_size) additional memory is necessary to build the candidates, which can be costly if the chromosome_size is very large.\n\nwarning: Warning\nAs it is in evolve!(), the decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing).\n\nwarning: Warning\nTHIS IS AN INTERNAL METHOD AND IT IS NOT MEANT TO BE USED DIRECTLY. IT IS CALLED FROM THE path_relink!() FUNCTION. Due to this reason, this method DOES NOT perform health checks on the arguments.\n\nArguments\n\nbrkga_data::BrkgaData: the BRKGA data.\nchromosome1::Array{Float64, 1} and chromosome2::Array{Float64, 1}: the chromosomes to be used to build the path.\naffect_solution::Function: function that takes two partial chromosomes / block of genes block1 and block2 and checks whether changing the keys from block1 to block2 affects the solution. For instance, suppose that the alleles/keys are used as threshold such that values > 0.5 activate a feature. Suppose we have block1 = [0.3, 0.4, 0.1] and block2 = [0.4, 0.1, 0.2]. Since all values are below 0.5, changing the keys from block1 to block2 does not chage the solution, and therefore, we can drop such change (and subsequentely decoding). The blocks can hold only one key/allele, sequential key blocks, of even the whole chromosome. affect_solution takes two views/subarrays. The function MUST HAVE the following signature\naffect_solution(block1::SubArray{Float64, 1},\n                block2::SubArray{Float64, 1})::Bool\nnote: Note\nThis function depends on the problem structure and how the   keys/alleles are used.\nblock_size::Int64: (posite) number of alleles to be exchanged at once in each iteration. If block_size == 1, the traditional path relinking is performed.\nmax_time::Float64: abort path-relinking when reach max_time. If max_time <= 0, no limit is imposed. Given in seconds.\npercentage::Float64: define the size, in percentage, of the path to build. Range [0, 1].\n\nReturns\n\nArray{Any, 1}: the best pair [fitness, chromosome] found during the relinking. If the relink is not possible due to homogeneity, -Inf returns in case of maximization, and Inf in case of minimization.\n\n\n\n\n\n"
},

{
    "location": "api/#BrkgaMpIpr.permutation_based_path_relink!",
    "page": "Library",
    "title": "BrkgaMpIpr.permutation_based_path_relink!",
    "category": "function",
    "text": "function permutation_based_path_relink!(brkga_data::BrkgaData,\n                                        chromosome1::Array{Float64, 1},\n                                        chromosome2::Array{Float64, 1},\n                                        affect_solution::Function,\n                                        block_size::Int64,\n                                        max_time::Float64,\n                                        percentage::Float64\n    )::Tuple{Float64, Array{Float64, 1}}\n\nPerform the permutation-based path relinking. In this method, the permutation induced by the keys in the guide solution is used to change the order of the keys in the permutation induced by the base solution.\n\nThe API will call decode!() function always with writeback = false. The reason is that if the decoder rewrites the chromosome, the path between solutions is lost and inadvertent results may come up. Note that at the end of the path relinking, the method calls the decoder with writeback = true in the best chromosome found to guarantee that this chromosome is re-written to reflect the best solution found.\n\nThis method is a multi-thread implementation. Instead of to build and decode each chromosome one at a time, the method builds a list of candidates, altering the alleles/keys according to the guide solution, and then decode all candidates in parallel. Note that O(chromosome_size^2 / block_size) additional memory is necessary to build the candidates, which can be costly if the chromosome_size is very large.\n\nwarning: Warning\nAs it is in evolve!(), the decoding is done in parallel using threads, and the user must guarantee that the decoder is THREAD-SAFE. If such property cannot be held, we suggest using single thread by setting the environmental variable JULIA_NUM_THREADS = 1 (see Julia Parallel Computing).\n\nwarning: Warning\nTHIS IS AN INTERNAL METHOD AND IT IS NOT MEANT TO BE USED DIRECTLY. IT IS CALLED FROM THE path_relink!() FUNCTION. Due to this reason, this method DOES NOT perform health checks on the arguments.\n\nArguments\n\nbrkga_data::BrkgaData: the BRKGA data.\nchromosome1::Array{Float64, 1} and chromosome2::Array{Float64, 1}: the chromosomes to be used to build the path.\naffect_solution::Function: not used in this function but kept to API compatibility.\nblock_size::Int64: not used in this function but kept to API compatibility.\nmax_time::Float64: abort path-relinking when reach max_time. If max_time <= 0, no limit is imposed. Given in seconds.\npercentage::Float64: define the size, in percentage, of the path to      build. Range [0, 1].\n\nReturns\n\nArray{Any, 1}: the best pair [fitness, chromosome] found during the relinking. If the relink is not possible due to homogeneity, -Inf returns in case of maximization, and Inf in case of minimization.\n\n\n\n\n\n"
},

{
    "location": "api/#Major-helper-functions-1",
    "page": "Library",
    "title": "Major helper functions",
    "category": "section",
    "text": "evolve_population!\ndirect_path_relink!\npermutation_based_path_relink!"
},

{
    "location": "api/#Index-1",
    "page": "Library",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"api.md\"]"
},

{
    "location": "contributing/#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "We aim to write an efficient and consistent code. So, you should keep in mind the balance between memory utilization and code efficiency and pay special attention to cache utilization. This aspect is very important when using multi-threads applications with shared memory."
},

{
    "location": "contributing/#Style-1",
    "page": "Contributing",
    "title": "Style",
    "category": "section",
    "text": "Please, follow the general Julia coding style. Since it is too long to describe all details here, study the code already written. However, in general,Name classes, methods, and variables as clear and meaningful as possible;\nWrite short commentaries on the code flow to reading more accessible and faster;\nProperly document the code, especially the data structures and methods definitions. Do not forget to link/refer them;\n4-space indentation, no trailing spaces, no tabs, Unix/POSIX end of line. Try to keep line within 80 columns and do not exceed 90 columns;\nDo not use one-liner branches. Always use if...end even it uses two more lines. The code must be as clear and easy to read as possible: # Don\'t do it\na > 1 && b += func(a)\n\n# Ah, way better and clear\nif a > 1\n    b += func(a)\nendAvoid dense expressions where possible e.g. prefer nested ifs over complex nested ?s;\nExplicit return should be preferred except in short-form method definitions. For functions that do not return values, explicitly use nothing at the end;\nDo not use system specific code/headers. Your code must compile in several systems with minimum change;\nMake sure that, for each method/function, you write unit tests that cover all corner cases, and few regular cases (> 1);\nDo not commit or do a pull request until the code pass in all tests."
},

]}
