decision-support-cvrp
=====================

Custom implementation of heuristics and meta-heuristics for the Capacitated Vehicle Routing Problem (CVRP): e.g. Clarke &amp; Wright Savings algorithm (CWS), Simulated Annealing (SA) and Genetic Algorithms (GAs).

# Building

Simply run `make` in the root dir.

# Basic usage

To run the CWS algorithm (parallel version) over dataset `doc/instances/A-n32-k5.vrp` run the following command:
`$ bin/cvrp --dataset=doc/instances/A-n32-k5.vrp --cws_version p`

You should get the following output:

    ROUTE LIST:

        ROUTE # 1:   0| 26|  6| 18| 28|  4| 11|  8|  9| 22| 14|  0,  [D:265][QS: 92][NS: 10]

        ROUTE # 2:   0| 16|  7| 13|  1| 12|  0,  [D:114][QS: 90][NS:  5]

        ROUTE # 3:   0| 20|  5| 25| 10| 15| 29| 27|  0,  [D:200][QS: 91][NS:  7]

        ROUTE # 4:   0| 23|  2|  3| 17| 19| 31| 21|  0,  [D:218][QS: 99][NS:  7]

        ROUTE # 5:   0| 24| 30|  0,  [D: 68][QS: 38][NS:  2]


    --------------------------------------------------------------------------------
    DATASET : doc/instances/A-n32-k5.vrp TOTAL COST (CWS) :  865, EXEC. TIME : 0.186
    --------------------------------------------------------------------------------


    SOLUTION VALIDATION REPORT: 
        [NUM STATIONS] = 32 32
        [STATION INCLUSION] = OK
        [MEAS. TOTAL DISTANCE] = 865 (SAVED IS 865)

    ROUTE LIST:

        ROUTE # 1:   0| 26|  6| 18| 28|  4| 11|  8|  9| 22| 14|  0,  [D:265][QS: 92][NS: 10]

        ROUTE # 2:   0| 16|  7| 13|  1| 12|  0,  [D:114][QS: 90][NS:  5]

        ROUTE # 3:   0| 20|  5| 25| 10| 15| 29| 27|  0,  [D:200][QS: 91][NS:  7]

        ROUTE # 4:   0| 23|  3|  2| 17| 19| 31| 21|  0,  [D:216][QS: 99][NS:  7]

        ROUTE # 5:   0| 24| 30|  0,  [D: 68][QS: 38][NS:  2]


    --------------------------------------------------------------------------------
    DATASET : doc/instances/A-n32-k5.vrp TOTAL COST (CWS) :  863, EXEC. TIME : 0.008
    --------------------------------------------------------------------------------


    SOLUTION VALIDATION REPORT: 
        [NUM STATIONS] = 32 32
        [STATION INCLUSION] = OK
        [MEAS. TOTAL DISTANCE] = 863 (SAVED IS 863)


# Available options

**-h, --help:** Here's the help page.

**--dataset <value>:** CVRP dataset to be fed to the algorithms. Should be in the same format as in http://www.coin-or.org/SYMPHONY/branchandcut/VRP/index.htm. e.g. --dataset=instances/A-n32-k5.vrp

**--cws_version <value>:** Version of the CWS algorithm to run (parallel or sequential). e.g. --cws_version=p. Use 'p' for running the parallel version, 's' for running the sequential version.

**--sa:** Apply a Simulated Annealing (SA) run on the solution provided by the CWS algorithm. By default SA is applied over the solution previously derived by the CWS algorithm. SA can work on a randomly generated solution if (1) the the --cws_version isn't used together with --sa or (2) if the --sa-random option is used.

**--sa_random:** Bypass any solution previously derived by the CWS algorithm and run SA over a randomly generated solution.

**--sa_cooling <value>:** Cooling ratio parameter for SA. Default is 0.005

**--sa_init_temp <value>:** Initial value for the temperature parameter for SA. If omitted, the method used in Osman1993 is applied to determine it.

**--sa_termination <value>:** Termination condition for SA, i.e. the number of neighbor search cycles performed by the SA method for each temperature epoch. Default is 25000.

**--sa_halting <value>:** Halting criterion parameter for SA, i.e. the number of max. temperature epochs a current solution determined via SA is allowed to remain unchanged. The breach of this condition effectively stops the SA procedure. Default is 300.

**--ga:** Run a Genetic Algorithm (GA) over a population generated after the results given by the CWS algorithm. By default GA considers in its initial population the solutions generated by the CWS algorithm (parallel or sequential, according to the --cws_version option) and CWS + SA (if --sa is used). To bypass this and create an initial population completely formed by random individuals, use the --ga_random option.

**--ga_random:** Bypass any solution previously derived by the CWS or SA algorithm and run GA over a population of randomly generated solutions.

**--ga_generation <value>:** Number of generations to be considered by the GA. Default is 10000.

**--ga_population <value>**: Population size to be considered by the GA. Default is 120.

**--ga_mutation <value>:** Mutation probability to be considered by the GA. Use values between 0.0 and 1.0 (inclusive). Default is 0.5.

**--matlab:** Create a .m file for graphical display of the final solution on MATLAB.
