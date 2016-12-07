#include "SA.h"

//! SA constructor.
/*!
*/
SA::SA()
{
}

//! SA destructor.
/*!
*/
SA::~SA()
{
}

//! Generates a random number (double floating point precision) in the interval 
//! [0.00, 1.00[.
/*!
*/
double Rand() {

    double r = (double) rand() / (double) (RAND_MAX + 1.0);

    return r;
}

double CalcAcceptanceProb(
        long delta, 
        double temp) {

    // If the route distance from the candidate solution is better, return 
    // a probability of 1.0, i.e. force the SA algorithm to accept this new 
    // solution.
    if (delta >= 0) {

        return 1.0;
    }

    // Calculate the probability of acceptance, using formula (6.1) on 
    // page 130 of Toth2002.
    double p = exp(((double) (delta)) / temp);

    return p;
}

//! Initialization of the cooling schedule parameters, 
//! as defined in Osman1993, i.e. determine an initial temperature value by 
//! exploring (without performing changes) the entire 1-interchange 
//! neighborhood of S_c.
/*!
*/
SA::CoolingScheduleParams CalcCoolingSchedule(CVRP * graph) {

    // Starting and final temperatures.
    uint32_t T_s = 0, T_f = UINT32_MAX;

    // Total number of feasible exchanges.
    uint32_t nfeas = 0;

    // Auxiliary variable for keeping track of changes.
    uint32_t abs_change = 0;

    // Other auxiliary variables;
    long dist_i_j   = 0;
    long dist_pj_j  = 0;
    long dist_pi_i  = 0;
    long dist_pi_j  = 0;
    long dist_pj_nj = 0;
    long dist_j_nj  = 0;
    long dist_pj_i  = 0;
    long dist_pi_ni = 0;
    long dist_i_ni  = 0;
    //long dist_j_i   = 0;

    long op_ij = 0;
    long op_ji = 0;
    long op_swap = 0;

    int qty_supplied_i = 0;
    int qty_supplied_j = 0;

    std::list <CVRP::Route *>::iterator it;
    std::list <CVRP::Route *>::iterator jt;

    // Go through all route combinations.
    for (it = graph->routes.begin(); it != graph->routes.end(); it++) {

        // For clearer code and readability, let's operate over an array, in 
        // which we could work with indexes.
        CVRP::Station** stations_i = graph->CreateRouteRep(*it);

        jt = it;
        jt++;

        for ( ; jt != graph->routes.end(); jt++) {

            CVRP::Station** stations_j = graph->CreateRouteRep(*jt);

            // Evaluate the gains from swapping 1 stations between 
            // routes route_i and route_j. Fortunately, this is never larger 
            // than OSMAN_ICHANGE_MAX_LAMBDA.
            for (int i = 1; i < (*it)->num_stations; i++) {
                for (int j = 1; j < (*jt)->num_stations; j++) {

                    dist_pj_i  = CVRP::Euc2DistCalc(stations_j[j - 1], stations_i[i]);
                    dist_i_j   = CVRP::Euc2DistCalc(stations_i[i], stations_j[j]);
                    dist_pj_j  = CVRP::Euc2DistCalc(stations_j[j - 1], stations_j[j]);

                    dist_pi_ni = CVRP::Euc2DistCalc(stations_i[i - 1], stations_i[i + 1]);
                    dist_pi_i  = CVRP::Euc2DistCalc(stations_i[i - 1], stations_i[i]);
                    dist_i_ni  = CVRP::Euc2DistCalc(stations_i[i], stations_i[i + 1]);

                    dist_pi_j  = CVRP::Euc2DistCalc(stations_i[i - 1], stations_j[j]);
                    //dist_j_i   = dist_i_j;

                    dist_pj_nj = CVRP::Euc2DistCalc(stations_j[j - 1], stations_j[j + 1]);
                    dist_j_nj  = CVRP::Euc2DistCalc(stations_j[j], stations_j[j + 1]);

                    // Check if the exchange is feasible, i.e. if the 
                    // route's qty_supplied value does not exceed the 
                    // capacity of the vehicles.
                    qty_supplied_j = 
                        (*jt)->qty_supplied 
                        + (stations_i[i]->demand);

                    // Change in distance caused by including point i in 
                    // route_j, i.e. testing operator (1,0)
                    if ((qty_supplied_j > graph->capacity)) {

                        op_ij = 0;

                    } else {

                        // A feasible (potential) change.
                        nfeas++;

                        op_ij = 
                            // Insert station i in between j - 1 and j.
                            + dist_pj_i
                            + dist_i_j
                            - dist_pj_j

                            // New distance between stations (i - 1) and (i + 1).
                            + dist_pi_ni
                            - dist_pi_i
                            - dist_i_ni;

                        abs_change = abs(op_ij);

                        if (abs_change > T_s)
                            T_s = abs_change;

                        if (abs_change < T_f)
                            T_s = abs_change;
                    }

                    qty_supplied_i = 
                        (*it)->qty_supplied 
                        + (stations_j[j]->demand);

                    if (qty_supplied_i > graph->capacity) {

                        op_ji = 0;

                    } else {

                        nfeas++;

                        // Change in distance caused by including point j in 
                        // route_i, i.e. testing operator (0,1).
                        op_ji = 
                            // Insert station j in between i - 1 and i.
                            + dist_pi_j
                            + dist_i_j
                            - dist_pi_i

                            // New distance between stations (j - 1) and (j + 1).
                            + dist_pj_nj
                            - dist_pj_j
                            - dist_j_nj;

                        abs_change = abs(op_ji);

                        if (abs_change > T_s)
                            T_s = abs_change;

                        if (abs_change < T_f)
                            T_s = abs_change;
                    }

                    if ((qty_supplied_i - stations_i[i]->demand) > graph->capacity
                        || (qty_supplied_j - stations_j[j]->demand) > graph->capacity) {

                        op_swap = 0;

                    } else {

                        // Update the number of feasible exchanges by 2 (i.e. 
                        // exchanging points i and j between routes i and j).
                        nfeas = nfeas + 2;

                        // Change in distance by exchanging points i and j between 
                        // route_i and route_j, i.e. testing operator (1,1).
                        op_swap = 
                            // Change in distance caused by including point i in 
                            // route_j, between stations j - 1 and j.
                            + dist_pj_i
                            + CVRP::Euc2DistCalc(stations_i[i], stations_j[j + 1])
                            - dist_pi_i
                            - dist_i_ni

                            // Change in distance caused by including point j in 
                            // route_i, between stations i - 1 and i.
                            + dist_pi_j
                            + CVRP::Euc2DistCalc(stations_j[j], stations_i[i + 1])
                            - dist_pj_j
                            - dist_j_nj;

                        abs_change = abs(op_swap);

                        if (abs_change > T_s)
                            T_s = abs_change;

                        if (abs_change < T_f)
                            T_f = abs_change;
                    }
                }
            }

            free(stations_j);
        }

        free(stations_i);
    }

    // Fill the cooling schedule parameters structure.
    SA::CoolingScheduleParams params;

    params.T_s = (double) T_s;
    params.T_f = (double) T_f;
    params.alpha = (graph->num_stations * nfeas);
    params.gamma = graph->num_stations;
    params.R = 3;

    printf("OSMAN'S 1-INTERCHANGE COOLING SCHEDULE METHOD: "\
            "Cooling Schedule Parameters: "\
            "\n\t[TS] = %0.1f"\
            "\n\t[TF] = %0.1f"\
            "\n\t[ALHPHA] = %d"\
            "\n\t[GAMMA] = %d"\
            "\n\t[R] = %d\n\n", 
        params.T_s,
        params.T_f,
        params.alpha,
        params.gamma,
        params.R);

    return params;
}

SA::SAPoint * SaveSAPoint(
        uint32_t k,
        uint32_t i,
        double T_k,
        long S_c,
        long S_n,
        long S_n_,
        long S_b,
        std::list <SA::SAPoint *> & sa_points) {

    SA::SAPoint * sa_point = new SA::SAPoint();

    sa_point->k = k;
    sa_point->i = i;
    sa_point->T_k = T_k;
    sa_point->S_c = S_c;
    sa_point->S_n = S_n;
    sa_point->S_n_ = S_n_;
    sa_point->S_b = S_b;

    sa_points.push_back(sa_point);

    return sa_point;
}

//! Prints a summary of a Evolutionary algorithm run, in the same fashion as 
//! for the CWS algorithm.
/*!
    Variables included in the summary are:
        - DATA: Dataset file the SA algorithm works on.
        - T_i: Initial temperature calculated via CalcCoolingSchedule().
        - A: ALPHA parameter or cooling rate.
        - S_i: Evaluation function value (i.e. total distance) of the initial 
                solution used by the SA algorithm.
        - S_b: Evaluation function value (i.e. total distance) of the best 
                solution found by the SA algorithm.
        - RT: SA algorithm running time (in milliseconds).
        - I: Total number of SA algorithm iterations.
*/
void PrintSAReport(
        const char * dataset, 
        SA::SAPoint * sa_point, 
        long S_i,
        double T_i,
        double alpha,
        uint32_t n, 
        double sa_time) {

    char report[512];

    sprintf(
            report, 
            "DATA: %s T_i: %0.4f A: %0.4f S_i: %ld S_b: %4ld RT: %0.3f I: %4d\n",
            dataset, 
            T_i,
            alpha,
            S_i,
            sa_point->S_b, 
            sa_time,
            (sa_point->k * n) + sa_point->i);

    printf("\n\n");
    for (unsigned int i = 0; i < (strlen(report) - 1); i++)
        printf("-");
    printf("\n");

    printf("%s", report);

    for (unsigned int i = 0; i < (strlen(report) - 1); i++)
        printf("-");
    printf("\n\n");
}

std::string GenerateSAFilename(std::string filename, std::string suffix) {

    char aux[16];

    sprintf(aux, "_%s", suffix.c_str());

    filename += aux + std::string(".csv");
    filename.replace(filename.find(".vrp"), 4, "");
    std::replace(filename.begin(), filename.end(), '-', '_');

    return filename;
}

void SAPoint2CSV(
        CVRP * graph,
        std::list <SA::SAPoint *> sa_points, 
        uint32_t n) {

    char aux_buffer[16];
    ofstream csv_file;

    std::string csv_filename = GenerateSAFilename(graph->dataset, "sa");

    printf("SA FILENAME: %s\n\n",
            csv_filename.c_str());

    csv_file.open(csv_filename.c_str());

    // Generate column names.
    csv_file << "Iteration, Temperature, Current, Candidate, BCandidate, Best\n";

    std::list <SA::SAPoint *>::iterator sat;

    for (sat = sa_points.begin(); sat != sa_points.end(); sat++) {

        csv_file << ((*sat)->k * n) + (*sat)->i << ", ";

        memset(aux_buffer, 0, sizeof(aux_buffer));

        sprintf(
                aux_buffer, 
                "%0.4f", 
                (*sat)->T_k);

        csv_file << aux_buffer << ", ";

        csv_file << ((*sat)->S_c) << ", " << ((*sat)->S_n) << ", " << ((*sat)->S_n_) << ", " << ((*sat)->S_b) << "\n";
    }

    csv_file.close();
}

void SA::Run(
        CVRP * graph, 
        double cooling_rate,
        double T_i,
        uint32_t L,
        uint16_t stop_n,
        SA::Phase start_phase,
        bool matlab) {

    // Current solution S_c.
    CVRP * S_c = new CVRP();

    // Current solution S_c.
    CVRP * S_n = new CVRP();

    // Best solution S_b.
    CVRP * S_b = new CVRP();

    // Start counting the execution time of the algorithm
    struct timeval tm;
    gettimeofday(&tm, NULL);
    double start_time = (tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0);

    // Check if the starting phase of the Osman's SA algorithm has been chosen 
    // to PHASE_1, i.e. it involves the calculation of an initial solution 
    // using the CWS algorithm + 2-opt exchange procedure (check CWS.cpp for 
    // details).
    if (start_phase == SA::PHASE_1) {

        // TODO:

    } else {

        // Set the initial solution S from the graph given as argument.
        S_c->CopyFrom(graph);
    }

    // Set S_b <- S_c (from the current solution S_c).
    S_b->CopyFrom(S_c);
    S_n->CopyFrom(S_c);

    // Save the number of iterations for which S_c remains unchanged, i.e. it 
    // is accepted after the stochastic step.
    int S_c_unchanged = 0;

    // Although the Osman1993 method is not implemented here, its method for 
    // calculation of initial temperature was maintained.
    CoolingScheduleParams params = CalcCoolingSchedule(graph);

    // Variables for Cooling Schedule. FIXME: We're using the initial 
    // temperature calculation method from Osman1993.
    // Use the formula set in our paper for determining the initial 
    // temperature 
    double T_k = (params.T_s / abs(log(0.5))) * 10.0;
    double T_s = T_k;

    double randomness = 0.0;
    double prob = 0.0;

    // List of SA iteration points, in order to keep track of the algorithm's 
    // performance. It will be later used to print the typical SA charts on 
    // MATLAB or Excel
    std::list <SAPoint *> sa_points;
    uint32_t sa_point_interval = ((L / MAX_SA_POINT) + 1);

    uint32_t k = 0;
    uint32_t i = 0;

    long best_distance = INT_MAX;

    // Keep 1-Interchange information, without performing changes. 
    // Notice that changes are always evaluated and performed over S_c.
    CVRP::IChange ichange;

    do {

        for (i = 0; i < L; i++) {

            // Reset ichange
            ichange.i = 0;
            ichange.j = 0;
            ichange.distance = 0;
            ichange.operation = 0;

            // Perform a random 1-interchange on S_n.
            ichange = S_n->RandomOneInterchange();

            // Keep track of the best value of S_n ever recorded, even if it 
            // is worse than S_c or S_b.
            if (S_n->total_distance < best_distance) {

                best_distance = S_n->total_distance;

                SaveSAPoint(
                        k, 
                        i, 
                        T_k, 
                        S_c->total_distance, 
                        S_n->total_distance,
                        best_distance, 
                        S_b->total_distance,
                        sa_points);
            }

            /*printf("SA: At iteration i = %d: "\
                    "\n\t[T_k] = %0.4f"\
                    "\n\t[S_i] = %ld"\
                    "\n\t[S_c] = %ld"\
                    "\n\t[S_n] = %ld"\
                    "\n\t[ichange.distance] = %ld"\
                    "\n\t[S_b] = %ld"\
                    "\n\t[S_c UNCHANGED FOR] = %d\n\n", 
                    i,
                    T_k,
                    graph->total_distance,
                    S_c->total_distance,
                    S_n->total_distance,
                    ichange.distance,
                    S_b->total_distance,
                    S_c_unchanged);*/

            // Save the info of this iteration (k,i). Only save aprox. 
            // MAX_SA_POINT values per iteration i.
            if ((i % sa_point_interval) == 0) {

                SaveSAPoint(
                        k, 
                        i, 
                        T_k, 
                        S_c->total_distance, 
                        S_n->total_distance,
                        best_distance, 
                        S_b->total_distance,
                        sa_points);
            }

            if (S_n->total_distance <= S_c->total_distance) {

                // Perform the changes on S_c.
                S_c->OneInterchange(ichange);

                S_c_unchanged = 0;

                /*printf("SA: S_n < S_c at iteration i = %d: "\
                        "\n\t[T_k] = %0.4f"\
                        "\n\t[S_i] = %ld"\
                        "\n\t[S_c] = %ld"\
                        "\n\t[S_n] = %ld"\
                        "\n\t[ichange.distance] = %ld"\
                        "\n\t[S_b] = %ld"\
                        "\n\t[S_c UNCHANGED FOR] = %d\n\n", 
                        i,
                        T_k,
                        graph->total_distance,
                        S_c->total_distance,
                        S_n->total_distance,
                        ichange.distance,
                        S_b->total_distance,
                        S_c_unchanged);*/

                if (S_n->total_distance <= S_b->total_distance) {

                    /*printf("SA: S_n < S_b at iteration i = %d: "\
                            "\n\t[T_k] = %0.4f"\
                            "\n\t[S_i] = %ld"\
                            "\n\t[S_c] = %ld"\
                            "\n\t[S_n] = %ld"\
                            "\n\t[ichange.distance] = %ld"\
                            "\n\t[S_b] = %ld"\
                            "\n\t[S_c UNCHANGED FOR] = %d\n\n", 
                            i,
                            T_k,
                            graph->total_distance,
                            S_c->total_distance,
                            S_n->total_distance,
                            ichange.distance,
                            S_b->total_distance,
                            S_c_unchanged);*/

                    S_b->CopyFrom(S_n);

                    SaveSAPoint(
                            k, 
                            i, 
                            T_k, 
                            S_c->total_distance, 
                            S_n->total_distance,
                            best_distance, 
                            S_b->total_distance,
                            sa_points);
                }

            } else {

                prob = CalcAcceptanceProb(
                            (S_c->total_distance - S_n->total_distance), 
                            T_k);

                randomness = Rand();

                if (prob > randomness) {

                    // Accept S_n.

                    // Perform the changes on S_c.
                    S_c->OneInterchange(ichange);

                    S_c_unchanged = 0;

                    /*printf("SA: Accepted S_n by chance at iteration i = %d: "\
                            "\n\t[T_k] = %0.4f"\
                            "\n\t[S_i] = %ld"\
                            "\n\t[S_c] = %ld"\
                            "\n\t[S_n] = %ld"\
                            "\n\t[ichange.distance] = %ld"\
                            "\n\t[S_b] = %ld"\
                            "\n\t[S_c UNCHANGED FOR] = %d\n\n", 
                            i,
                            T_k,
                            graph->total_distance,
                            S_c->total_distance,
                            S_n->total_distance,
                            ichange.distance,
                            S_b->total_distance,
                            S_c_unchanged);*/

                } else {

                    // Do not accept S_n and keep S_c.

                    // FIXME: This is computationally expensive and should be 
                    // avoided.
                    S_n->CopyFrom(S_c);

                    // FIXME: In addition to the copy made above, I figured 
                    // I have to do this, just because the order of the 
                    // Station vector is different in S_n and S_c after the 
                    // previous instruction.
                    S_c->CopyFrom(S_n);

                    /*printf("SA: Keeping S_c at iteration i = %d: "\
                            "\n\t[T_k] = %0.4f"\
                            "\n\t[S_i] = %ld"\
                            "\n\t[S_c] = %ld"\
                            "\n\t[S_n] = %ld"\
                            "\n\t[ichange.distance] = %ld"\
                            "\n\t[S_b] = %ld"\
                            "\n\t[S_c UNCHANGED FOR] = %d\n\n", 
                            i,
                            T_k,
                            graph->total_distance,
                            S_c->total_distance,
                            S_n->total_distance,
                            ichange.distance,
                            S_b->total_distance,
                            S_c_unchanged);*/

                    S_c_unchanged++;
                }
            }
        }

        printf("SA: At iteration k = %d: "\
                "\n\t[T_k] = %0.4f"\
                "\n\t[S_i] = %ld"\
                "\n\t[S_c] = %ld"\
                "\n\t[S_n] = %ld"\
                "\n\t[S_b] = %ld"\
                "\n\t[S_c UNCHANGED FOR] = %d"\
                "\n\t[BEST S_n SO FAR] = %ld"\
                "\n\t[SA POINTS]: %ld\n\n",
                k,
                T_k,
                graph->total_distance,
                S_c->total_distance,
                S_n->total_distance,
                S_b->total_distance,
                S_c_unchanged,
                best_distance,
                sa_points.size());

        // Cool it down...
        T_k = T_k * (1 - cooling_rate);

        // Increment L_k
        k++;

    } while (S_c_unchanged < stop_n);

    gettimeofday(&tm, NULL);
    double sa_time = ((tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0)) - start_time;

    // Save the last SAPoint.
    SAPoint * sa_point = SaveSAPoint(
            k, 
            i, 
            T_k, 
            S_c->total_distance, 
            S_n->total_distance,
            best_distance, 
            S_b->total_distance,
            sa_points);

    S_b->PrintRoutes();
    S_b->Validate();

    graph->CopyFrom(S_b);

    // Print a SA report.
    PrintSAReport(
            S_b->dataset.c_str(), 
            sa_point, 
            graph->total_distance, 
            T_s, 
            cooling_rate, 
            L, 
            sa_time);

    // Save a .csv file with the SA results (SAPoint).
    /*SAPoint2CSV(
        S_b,
        sa_points, 
        L);*/

    // Save a nice graph for the final SA result.
    Matlab::GenerateCVRPGraph(S_b, Matlab::DRAW_ROUTES, "sa");

    // Free the memory allocated by each one of the SAPoint instances.
    std::list <SAPoint *>::iterator sat;

    for (sat = sa_points.begin(); sat != sa_points.end(); ) {

        delete (*sat);
        sat = sa_points.erase(sat);
    }
}


