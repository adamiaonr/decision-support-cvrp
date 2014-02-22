#include "CWS.h"
#include "Quicksort.h"

//! CWS constructor.
/*!
*/
CWS::CWS()
{
}

//! CWS destructor.
/*!
*/
CWS::~CWS()
{
}

void PrintSavings(CWS::Saving * savings, int dimension) {

    int s_dist = 0, s_sav = 0;

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {

            if (j > i) {

                printf(
                        "%4ld|", 
                        savings[s_sav++].s_ij);

            } else if (j < i) {

                printf(
                        "%4ld|", 
                        savings[s_dist++].dist_ij);

            } else {

                printf("    |");

            }
        }

        printf("\n");
    }
}

int UpdateRoute(
        CVRP * graph,
        CVRP::Route * dst,
        CVRP::Route * src,
        CWS::MergePosition merge_pos,
        CWS::Saving s) {

    // This will require a series of updates...

    // Total demand satisfied along the route.
    dst->qty_supplied += src->qty_supplied;

    // Distance, taking the saving value between (i,j) into account.
    dst->distance += src->distance - s.s_ij;

    // Total graph's distance update.
    graph->total_distance += - s.s_ij;

    // Number of stations.
    dst->num_stations += src->num_stations;

    // Start and end stations of the route.
    switch (merge_pos) {

        case CWS::LEFT:

            dst->start_station = src->start_station;
            break;

        case CWS::RIGHT:

            dst->end_station = src->end_station;
            break;

        default:

            // An error occurred, return -1.
            return -1;
    }

    // Update the RouteInfo pointers on the stations of route src.
    CVRP::Station * curr_station = src->end_station;
    //CVRP::Station * depot = graph->GetDepot();

    while (curr_station->type != CVRP::DEPOT) {

        curr_station->route = dst;
        curr_station = curr_station->prev;
    }

    // Delete the src Route and remove it from the list of routes.
    graph->DeleteRoute(src);

    // Return the dst route index.
    return dst->index;
}

int MergeRoutes(
        CVRP * graph,
        CWS::Saving & s,
        int capacity) {

    CVRP::Station * i = s.station_i;
    CVRP::Station * j = s.station_j;

    //int route_i_index = i->route->index;
    //int route_j_index = j->route->index;

    // Are the stations in the same route? If so, abort.
    if (i->route->index == j->route->index) {

        /*printf("MERGE: Stations %3d and %3d in the same route %2d. ABORTING.\n", 
                i->index, j->index, i->route->index);*/

        return -1;
    }

    int route_i_qty_supplied = i->route->qty_supplied;
    int route_j_qty_supplied = j->route->qty_supplied;

    // Is the capacity constraint violated by this merge operation?
    if ((route_i_qty_supplied + route_j_qty_supplied) > capacity) {

        /*printf("MERGE: Addition of link (%3d, %3d) violates capacity constraint "\
                "for routes %2d or %2d (%3d + %3d = %3d > %3d). ABORTING.\n", 
                i->index, j->index, 
                i->route->index, j->route->index,
                i->route->qty_supplied, j->route->qty_supplied, 
                i->route->qty_supplied + j->route->qty_supplied, 
                (int) capacity);*/

        return -1;
    }

    // Get a vertex descriptor for the depot.
    //CVRP::Station * depot = graph->GetDepot();

    int route_i_size = i->route->num_stations;
    int route_j_size = j->route->num_stations;
    int route_index = -1;

    if (graph->GetNextStation(i)->type == CVRP::DEPOT 
            && graph->GetPrevStation(j)->type == CVRP::DEPOT) {

        // This will make a route disappear. Always make the smaller 
        // route disappear. If equal, route containing j disappears.
        if (route_i_size < route_j_size) {

            // Update the RouteInfo object associated with j, according to 
            // the info on i's RouteInfo (this includes the elimination of 
            // i's RouteInfo).
            route_index = UpdateRoute(
                    graph,
                    j->route, 
                    i->route, 
                    CWS::LEFT,
                    s);

        } else {

            route_index = UpdateRoute(
                    graph,
                    i->route, 
                    j->route, 
                    CWS::RIGHT,
                    s);
        }

        // Simply update the graph, adding a new Edge (i,j) and deleting 
        // the previous (0,j) and (i,0).
        graph->LinkStations(i, j);

    } else if (graph->GetNextStation(j)->type == CVRP::DEPOT 
            && graph->GetPrevStation(i)->type == CVRP::DEPOT) {

        if (route_i_size < route_j_size) {

            // Update the RouteInfo object associated with j, according to 
            // the info on i's RouteInfo (this includes the elimination of 
            // i's RouteInfo).
            route_index = UpdateRoute(
                    graph,
                    j->route, 
                    i->route, 
                    CWS::RIGHT,
                    s);

        } else {

            route_index = UpdateRoute(
                    graph,
                    i->route, 
                    j->route, 
                    CWS::LEFT,
                    s);
        }

        // Simply update the graph, adding a new Edge (j,i) and deleting 
        // the previous (0,i) and (j,0).
        graph->LinkStations(j, i);

    } else {

        /*printf("MERGE: At least one of the nodes %3d and %3d is interior to a "\
                "route. ABORTING.\n", i->index, j->index);*/
    }

    return route_index;
}

int SeqMergeRoutes(
        CVRP * graph,
        CVRP::Route * route,
        CWS::Saving & s,
        int capacity) {

    CVRP::Station * linking_station = NULL;
    CVRP::Station * candidate_station = NULL;

    if (!s.active) {

        return -1;
    }

    if (s.station_i->type == CVRP::DEPOT 
            || s.station_j->type == CVRP::DEPOT) {

        return -1;
    }

    if (s.station_i->route->index == route->index) {

        linking_station = s.station_i;
        candidate_station = s.station_j;

    } else {

        if (s.station_j->route->index == route->index) {

            linking_station = s.station_j;
            candidate_station = s.station_i;

        } else {

            return -1;
        }
    }

    // Are the nodes in the same route?
    if (linking_station->route->index == candidate_station->route->index) {

        /*printf("SEQMERGE: Nodes %2d and %2d in the same route %2d. ABORTING.\n", 
                linking_station->index, 
                candidate_station->index, 
                linking_station->route->index);*/

        return -1;
    }

    // Is the capacity constraint violated by this merge?
    if ((route->qty_supplied + candidate_station->route->qty_supplied) > capacity) {

        /*printf("SEQMERGE: Addition of link (%3d, %3d) violates capacity constraint "\
                "for routes %2d or %2d (%3d + %3d = %3d > %3d). ABORTING.\n", 
                s.station_i->index, s.station_j->index, 
                s.station_i->route->index, s.station_j->route->index,
                s.station_i->route->qty_supplied, s.station_j->route->qty_supplied, 
                s.station_i->route->qty_supplied + s.station_j->route->qty_supplied, 
                (int) capacity);*/

        return -1;
    }

    int route_index = -1;
    CVRP::Station * depot = graph->GetDepot();

    if (graph->GetNextStation(linking_station) == depot 
            && graph->GetPrevStation(candidate_station) == depot) {

        /*printf("SEQMERGE: Condition (i,0) (0,j) satisfied. Merging route "\
                "%2d to %2d via link (%2d,%2d)\n", 
                candidate_station->route->index, linking_station->route->index, 
                linking_station->index, candidate_station->index);*/

        route_index = UpdateRoute(
                    graph,
                    linking_station->route, 
                    candidate_station->route, 
                    CWS::RIGHT,
                    s);

        graph->LinkStations(linking_station, candidate_station);

        s.active = false;

    } else if (graph->GetPrevStation(linking_station) == depot 
            && graph->GetNextStation(candidate_station) == depot) {

        /*printf("SEQMERGE: Condition (j,0) (0,i) satisfied. Merging route "\
                "%2d to %2d via link (%2d,%2d)\n", 
                candidate_station->route->index, linking_station->route->index, 
                candidate_station->index, linking_station->index);*/

        route_index = UpdateRoute(
                    graph,
                    linking_station->route, 
                    candidate_station->route, 
                    CWS::LEFT,
                    s);

        graph->LinkStations(candidate_station, linking_station);

        s.active = false;

    } else {

        /*printf("SEQMERGE: At least one of the nodes %2d and %2d is interior to a "\
                "route. ABORTING.\n", 
                linking_station->index, 
                candidate_station->index);*/
    }

    return route_index;
}

/*int CWS::TwoOptExchange(CVRP * graph, CVRP::Route * route) {

    long best_distance = route->distance;
    long new_distance = INT_MAX;

    // For clearer code and readability, let's operate over an array, in 
    // which we could work with indexes.
    CVRP::Station** route_stations = (CVRP::Station**) calloc(route->num_stations + 2, sizeof(CVRP::Station*));
    CVRP::Station * curr_route_station = route->start_station;

    // Just a reference to the graph's depot.
    CVRP::Station * depot = graph->GetDepot();

    for (int i = 1; i < route->num_stations + 2; i++) {

                route_stations[i] = curr_route_station;
                curr_route_station = curr_route_station->next;
    }

    // Set the edges of the array to point to the depot.
    route_stations[0] = depot;
    route_stations[route->num_stations + 1] = depot;

    CVRP::Station * aux = NULL;

    // Start after the depot, hence i = 1.
    for (int i = 1; i < route->num_stations - 1; i++) {
        for (int j = (i + 1); j < route->num_stations; j++) {

            // Note that we're not using the next or prev pointers here, 
            // one is just evaluating what would happen if 2 stations 
            // swapped.
            new_distance = 
                    + route->distance 
                    + CVRP::Euc2DistCalc(route_stations[i - 1], route_stations[j])
                    + CVRP::Euc2DistCalc(route_stations[i], route_stations[j + 1])
                    - CVRP::Euc2DistCalc(route_stations[i - 1], route_stations[i])
                    - CVRP::Euc2DistCalc(route_stations[j], route_stations[j + 1]);

            if (new_distance < best_distance) {

                // Swap the stations in the graph.
                graph->SwapStations(route_stations[i], route_stations[j], j - i);

                // Swap the stations in the local array representing the 
                // Route.
                int m = i;

                for (int k = j; k > m; k--) {

                    aux = route_stations[k];
                    route_stations[k] = route_stations[m];
                    route_stations[m++] = aux;
                }

                best_distance = new_distance;
            }
        }
    }

    free(route_stations);

    return best_distance;
}*/

int CalcSavings(CVRP * graph, CWS::Saving * savings) {

    long savings_index = 0;
    long dist_ij = 0;
    long dist_Di = 0;
    long dist_Dj = 0;

    CVRP::Station * depot = graph->GetDepot();

    std::vector <CVRP::Station *>::iterator it;
    std::vector <CVRP::Station *>::iterator jt;

    it = graph->stations.begin();
    it++;

    for ( ; it != graph->stations.end(); it++) {

        // Distance between depot and station i.
        dist_Di = CVRP::Euc2DistCalc(
                depot->coords, 
                (*it)->coords);

        jt = it;
        jt++;

        for ( ; jt != graph->stations.end(); jt++) {

            /*printf("\t[INDEXES]: (%3d, %3d)\n", 
                    (*it)->index, 
                    (*jt)->index);*/

            // Distance between stations i and j.
            dist_ij = CVRP::Euc2DistCalc(
                    (*it)->coords, 
                    (*jt)->coords);

            // Distance between depot and station j.
            dist_Dj = CVRP::Euc2DistCalc(
                    depot->coords, 
                    (*jt)->coords);

            // Calculate the savings s(i,j), save it on the savings array, 
            // passed as argument.
            savings[savings_index].i = (*it)->index - 1;
            savings[savings_index].station_i = *it;
            savings[savings_index].j = (*jt)->index - 1;
            savings[savings_index].station_j = *jt;
            savings[savings_index].dist_ij = dist_ij;
            savings[savings_index].s_ij = (dist_Di + dist_Dj - dist_ij);
            savings[savings_index].active = true;

            /*printf("%4ld : (%3d,%3d) : %4ld (%3d)\n", 
                    savings_index,
                    savings[savings_index].i, 
                    savings[savings_index].j, 
                    savings[savings_index].s_ij,
                    (int) savings[savings_index].dist_ij);*/

            savings_index++;
        }
    }

    return 0;
}

void InitRoutes(CVRP * graph) {

    CVRP::Station * depot = graph->GetDepot();

    std::vector <CVRP::Station *>::iterator it;

    it = graph->stations.begin();
    it++;

    for ( ; it != graph->stations.end(); it++) {

        // Create the initial set of routes (0,i,0), i.e. add two edges 
        // (depot, i) and (i, depot).
        graph->LinkStations(depot, *it);
        graph->LinkStations(*it, depot);

        // For each new route (0,i,0) add a new Route instance to graph.
        CVRP::Route * route_info = graph->AddRoute(
                (*it)->demand,
                2 * CVRP::Euc2DistCalc(
                        depot->coords, 
                        (*it)->coords),
                (*it),
                (*it));

        // Assign the Route instance pointed by route_info to the current 
        // Station i. This way, when merging routes via end Stations i and j, 
        // one can quickly access information of the routes to which i and j 
        // belong.
        graph->Assign2Route(*it, route_info);
    }
}

void PrintCWSReport(const char * dataset, int total_cost, double cws_time) {

    char report[256];

    sprintf(
            report, 
            "DATASET : %s TOTAL COST (CWS) : %4d, EXEC. TIME : %0.3f\n",
            dataset, total_cost, cws_time);

    printf("\n\n");
    for (int i = 0; i < (strlen(report) - 1); i++)
        printf("-");
    printf("\n");

    printf("%s", report);

    for (int i = 0; i < (strlen(report) - 1); i++)
        printf("-");
    printf("\n\n");
}

void CWS::Run(
        CVRP * graph,
        CWS::Version version,
        bool matlab) {

    //graph->PrintStations();

    // Start counting the execution time of the algorithm
    struct timeval tm;
    gettimeofday(&tm, NULL);
    double start_time = (tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0);

    int number_nodes = graph->num_stations;

    long savings_size = ((number_nodes * number_nodes) - ((3 * number_nodes) - 2))/2;
    CWS::Saving * savings = (CWS::Saving *) calloc (savings_size, sizeof(CWS::Saving));

    // 1) Calculate the savings s(i, j) = d(D, i) + d(D, j) - d(i, j) for 
    // every pair (i, j) of demand points.
    CalcSavings(graph, savings);

    // 2) Print the savings list.
    //PrintSavings(savings, number_nodes);

    // 3) Sort the saving values in the savings list.
    Quicksort::Sort(savings, 0, savings_size - 1);

    // 4) Build an initial set of Routes, with the form 
    // (depot, Station i, depot).
    InitRoutes(graph);
    //graph->PrintRoutes();

    // 5)

    // 5.1) Parallel version of the Clarke and Wright Savings.

    int route_index = 0;

    if (version == CWS::PARALLEL) {

        for (long i = 0; i < (savings_size); i++) {

            route_index = MergeRoutes(
                                graph,
                                savings[i], 
                                graph->capacity);
        }

    } else {

        // 5.2) Sequential version of the Clarke and Wright Savings.

        std::list <CVRP::Route *>::iterator it;

        for (it = graph->routes.begin(); it != graph->routes.end(); it++) {

            for (long i = 0; i < (savings_size); i++) {

                if ((route_index = SeqMergeRoutes(
                    graph,
                    *it,
                    savings[i], 
                    graph->capacity)) >= 0) {

                    // Route extended, restart the savings list cycle.
                    i = 0;

                    savings[i].active = false;

                } else {

                    //CWS::print_route(&routes[route_index]);
                }
            }
        }
    }

    gettimeofday(&tm, NULL);
    double cws_time = ((tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0)) - start_time;

    int total_cost = 0;

    graph->PrintRoutes();

    std::list <CVRP::Route *>::iterator it;

    for (it = graph->routes.begin(); it != graph->routes.end(); it++) {

        total_cost += (*it)->distance;
    }

    PrintCWSReport(graph->dataset.c_str(), total_cost, cws_time);
    graph->Validate();

    Matlab::GenerateCVRPGraph(graph, Matlab::DRAW_ROUTES, "cws");

    gettimeofday(&tm, NULL);
    cws_time = ((tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0)) - start_time;

    total_cost = 0;

    for (it = graph->routes.begin(); it != graph->routes.end(); it++) {

        total_cost += graph->TwoOptExchange(*it, CVRP::HARD_2_OPT);
    }

    gettimeofday(&tm, NULL);
    double two_opt_time = ((tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0)) - start_time;

    graph->PrintRoutes();

    PrintCWSReport(
            graph->dataset.c_str(), 
            total_cost, 
            (two_opt_time - cws_time));

    graph->Validate();

    Matlab::GenerateCVRPGraph(graph, Matlab::DRAW_ROUTES, "2opt");

    free(savings);
}

