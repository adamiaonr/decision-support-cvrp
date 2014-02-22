#include "SA.h"
#include "CWS.h"
#include "CVRP.h"
#include "Evolutionary.h"
#include "argvparser.h"

using namespace std;
using namespace CommandLineProcessing;

ArgvParser * CreateArgvParser() {

    ArgvParser * cmds = new ArgvParser();

    /*cmds->setIntroductoryDescription("Clarke Wright Savings (CWS) algorithm "\
            "implementation for the CVRP by adamiaonr@gmail.com");*/
    cmds->setIntroductoryDescription("Custom implementation of heuristics and "\
            "meta-heuristics for the Capacitated Vehicle Routing Problem "\ 
            "(CVRP): e.g. Clarke & Wright Savings algorithm (CWS), Simulated "\
            "Annealing (SA) and Genetic Algorithms (GAs).\nImplementation by "\
            "adamiaonr@gmail.com");

    cmds->setHelpOption("h", "help", "Here's the help page.");

    cmds->defineOption(
            "dataset",
            "CVRP dataset to be fed to the algorithms. Should be in the same "\
                "format as in "\
                "http://www.coin-or.org/SYMPHONY/branchandcut/VRP/index.htm.\n"\
            "e.g. --dataset=instances/A-n32-k5.vrp",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "cws_version",
            "Version of the CWS algorithm to run (parallel or sequential).\n"\
            "e.g. for the parallel version use --cws_version=p\n"
            "Use 'p' for running the parallel version, 's' for running the "\
                "sequential version.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "sa",
            "Apply a Simulated Annealing (SA) run on the solution provided "\
                    "by the CWS algorithm. By default SA is applied over the "\
                    "solution previously derived by the CWS algorithm. SA "\
                    "can work on a randomly generated solution if (1) the "\
                    "the --cws_version isn't used together with --sa or (2) if "\
                    "the --sa-random option is used.",
            ArgvParser::NoOptionAttribute);

    cmds->defineOption(
            "sa_random",
            "Bypass any solution previously derived by the CWS algorithm and "\
                    "run SA over a randomly generated solution.",
            ArgvParser::NoOptionAttribute);

    cmds->defineOption(
            "sa_cooling",
            "Cooling ratio parameter for SA. Default is 0.005",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "sa_init_temp",
            "Initial value for the temperature parameter for SA. If "\
                    "omitted, the method used in Osman1993 is applied to "\
                    "determine it.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "sa_termination",
            "Termination condition for SA, i.e. the number "\
                    "of neighbor search cycles performed by the SA method "\
                    "for each temperature epoch. Default is 25000.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "sa_halting",
            "Halting criterion parameter for SA, i.e. the number "\
                    "of max. temperature epochs a current solution determined "\
                    "via SA is allowed to remain unchanged. The breach of this "\
                    "condition effectively stops the SA procedure. "\
                    "Default is 300.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "ga",
            "Run a Genetic Algorithm (GA) over a population generated after "\
                    "the results given by the CWS algorithm. By default GA "\
                    "considers in its initial population the solutions generated "\
                    "by the CWS algorithm (parallel or sequential, according to "\
                    "the --cws_version option) and CWS + SA (if --sa is used). "\
                    "To bypass this and create an initial population completely "\
                    "formed by random individuals, use the --ga_random option.",
            ArgvParser::NoOptionAttribute);

    cmds->defineOption(
            "ga_random",
            "Bypass any solution previously derived by the CWS or SA algorithm and "\
                    "run GA over a population of randomly generated solutions.",
            ArgvParser::NoOptionAttribute);

    cmds->defineOption(
            "ga_generation",
            "Number of generations to be considered by the GA. Default is 10000.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "ga_population",
            "Population size to be considered by the GA. Default is 120.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "ga_mutation",
            "Mutation probability to be considered by the GA. Use values "\
                    "between 0.0 and 1.0 (inclusive). Default is 0.5.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            "matlab",
            "Create a .m file for graphical display of the final "\
                    "solution on MATLAB.",
            ArgvParser::NoOptionAttribute);

    return cmds;
}

int main (int argc, char **argv) {

    // Initialize random number generator (required for Simulated Annealing 
    // (SA)).
    srand(time(NULL));

    // Boolean variables to handle the application's input options.
    bool is_cws_version = false;
    bool is_sa          = false;
    bool is_sa_random   = false;
    bool is_ga          = false;
    bool is_ga_random   = false;
    bool is_matlab      = false;

    // Set options to default values for SA.
    double sa_cooling       = SA_DEFAULT_COOLING_RATE;
    double sa_init_temp     = 512.0;
    uint32_t sa_termination = SA_DEFAULT_TERMINATION_CRITERION;
    uint16_t sa_halting     = SA_DEFAULT_HALTING_CRITERION;

    // Set options to default values for GA.
    uint32_t ga_generation  = GA_DEFAULT_GENERATIONS;
    long ga_population      = GA_DEFAULT_POPULATION_SIZE;
    double ga_mutation      = GA_DEFAULT_MUTATION_PROB;

    ArgvParser * cmds = CreateArgvParser();

    char * dataset = NULL;
    char * cws_version = NULL;
    CWS::Version _cws_version;

    int result = cmds->parse(argc, argv);

    if (result != ArgvParser::NoParserError) {

        fprintf(stderr, "%s\n", cmds->parseErrorDescription(result).c_str());

        if (result != ArgvParser::ParserHelpRequested) {
            fprintf(stderr, "Use the option -h for help.\n");
        }

        return -1;

    } else {

        if (cmds->foundOption("dataset")) {

            dataset = (char *) cmds->optionValue("dataset").c_str();

        } else {

            fprintf(stderr, "No dataset file specified. Use the option -h for help.\n");

            return -1;
        }

        if (cmds->foundOption("cws_version")) {

            cws_version = (char *) cmds->optionValue("cws_version").c_str();

            if (strcmp(cws_version, PARALLEL_VERSION) == 0) {

                _cws_version = CWS::PARALLEL;
                is_cws_version = true;

            } else if (strcmp(cws_version, SEQUENTIAL_VERSION) == 0) {

                _cws_version = CWS::SEQUENTIAL;
                is_cws_version = true;

            } else {

                fprintf(stderr, "Unrecognized CWS version option. Use the "\
                        "option -h for help.\n");

                return -1;
            }
        }

        if (cmds->foundOption("sa")) {

            is_sa = true;

            if (cmds->foundOption("sa_random")) {
                is_sa_random = true;
            }

            if (cmds->foundOption("sa_cooling")) {

                sa_cooling = 
                        (double) atof((char *) cmds->optionValue("sa_cooling").c_str());
            }

            if (cmds->foundOption("sa_init_temp")) {

                sa_init_temp = 
                        (double) atof((char *) cmds->optionValue("sa_init_temp").c_str());
            }

            if (cmds->foundOption("sa_termination")) {

                sa_termination = 
                        (uint32_t) atoi((char *) cmds->optionValue("sa_termination").c_str());
            }

            if (cmds->foundOption("sa_halting")) {
                sa_halting = 
                        (uint16_t) atoi((char *) cmds->optionValue("sa_halting").c_str());
            }
        }

        if (cmds->foundOption("ga")) {

            is_ga = true;

            if (cmds->foundOption("ga_random")) {
                is_ga_random = true;
            }

            if (cmds->foundOption("ga_generation")) {
                ga_generation = 
                        (uint32_t) atoi((char *) cmds->optionValue("ga_generation").c_str());
            }

            if (cmds->foundOption("ga_population")) {
                ga_population = 
                        (long) atoi((char *) cmds->optionValue("ga_population").c_str());
            }

            if (cmds->foundOption("ga_mutation")) {
                ga_mutation = 
                        (double) atof((char *) cmds->optionValue("ga_mutation").c_str());

                if (ga_mutation < 0.0 || ga_mutation > 1.0) {

                    fprintf(stderr, "Invalid 'ga_mutation' option value. Use the "\
                            "format '0.5' and keep it in the interval [0.0; 1.0].\n");

                    return -1;
                }
            }
        }

        if (cmds->foundOption("matlab")) {

            is_matlab = true;
        }
    }

    // Stop it from going ahead if the help flag was the one given as 
    // input.
    if (result == ArgvParser::ParserHelpRequested) {

        return -1;
    }

    // Pointers to different method's solutions.
    CVRP * cws_graph = NULL;
    CVRP * sa_graph = NULL;

    // Now, play along with the input options...

    // Run the _cws_version of the CWS algorithm.
    if (is_cws_version) {

        // Create an instance of a CVRP graph and fill with the info on the 
        // dataset file.
        cws_graph = new CVRP(std::string(dataset));
        CWS::Run(cws_graph, _cws_version, is_matlab);
    }

    // Run Simulated Annealing (SA)...
    if (is_sa) {

        // ...over a randomly generated solution or the previous solution 
        // determined via the CWS algorithm.
        if (is_sa_random || is_cws_version == false) {

            // FIXME: This is just a trick to determine the appropriate number 
            // of routes. It seem's quite lame, but as the 'optimal' number of 
            // routes isn't extracted from the dataset file, we have to do this 
            // for now (we use the CWS::SEQUENTIAL version as it always 
            // provided solutions with a number of routes equal to the 'optimal' 
            // value).
            CVRP * dummy = new CVRP(std::string(dataset));
            CWS::Run(dummy, CWS::SEQUENTIAL, false);

            if (CVRP::GenerateRandomIndividual(
                    sa_graph,
                    dummy->dataset,
                    (int) dummy->routes.size()) < 0) {

                fprintf(stderr, "Fatal error when generating random solution. "\
                        "Aborting.\n");

                return -1;
            }

        } else {

            sa_graph = new CVRP();
            sa_graph->CopyFrom(cws_graph);
        }

        SA::Run(
                sa_graph, 
                sa_cooling, 
                sa_init_temp, 
                sa_termination, 
                sa_halting, 
                SA::PHASE_2,
                is_matlab);
    }

    // Run the Genetic Algorithm (GA) implementation.
    if (is_ga) {

        // If a totally random population has been requested OR there are 
        // no CWS and SA solutions to work with, run GA over a completely 
        // random population.
        if (is_ga_random || (is_cws_version == false && is_sa == false)) {

            Evolutionary::Run(
                    std::string(dataset),
                    ga_generation,
                    ga_population,
                    ga_mutation,
                    is_matlab);

        } else {

            // Otherwise, add the CWS and SA solutions (one of them may 
            // be NULL)...
            Evolutionary::Run(
                    cws_graph,
                    sa_graph,
                    ga_generation,
                    ga_population,
                    ga_mutation,
                    is_matlab);
        }
    }

    return 1;
}

