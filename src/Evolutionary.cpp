#include "Evolutionary.h"

uint32_t global_id = 0;

//! Evolutionary constructor.
/*!
*/
Evolutionary::Evolutionary()
{
}

//! Evolutionary destructor.
/*!
*/
Evolutionary::~Evolutionary()
{
}

//! Prints a summary of a Evolutionary algorithm run, in the same fashion as 
//! for the CWS and SA algorithms.
/*!
*/
void PrintEvolutionaryReport(
        const char * dataset, 
        long S_i, 
        long S_b,
        uint32_t S_b_generation,
        int pop_size,
        double mutation_prob,
        uint32_t generations, 
        double ga_time) {

    char report[512];

    sprintf(
            report, 
            "DATA: %s S_i: %4ld S_b: %4ld @ %d P: %3d P(M): %0.4f G: %4d RT: %0.3f\n ",
            dataset, 
            S_i,
            S_b,
            S_b_generation,
            pop_size,
            mutation_prob, 
            generations,
            ga_time);

    printf("\n\n");
    for (int i = 0; i < (strlen(report) - 1); i++)
        printf("-");
    printf("\n");

    printf("%s", report);

    for (int i = 0; i < (strlen(report) - 1); i++)
        printf("-");
    printf("\n\n");
}

//! Sort the candidates vector in ascending order, i.e. lowest evaluation 
//! function values (i.e. total distance) first.
/*!
*/
bool SortCandidates(CVRP * a, CVRP * b) {

    return a->total_distance < b->total_distance;
}

void RandomTournamentSelection(
        std::vector <CVRP *> & contestants, 
        std::vector <CVRP *> & fittest) {

    long contestants_size = contestants.size();
    long pop_size = contestants_size;

    // 1) Sort the contestants by evaluation function value, in ascending 
    // order, i.e. smallest values first.
    sort(
            contestants.begin(), 
            contestants.end(), 
            SortCandidates);

    /*for (int m = 0; m < pop_size; m++) {

        printf("contestant[%d] = %ld\n", m, contestants[m]->total_distance);
    }*/

    // 2) Select the first, more valuable, pop_size / 10 contestants and make 
    // them automatically eligible to the next generation.
    for (int m = 0; m < (pop_size / 10); m++) {

        // 2.1) Send it to Noah's Ark...
        fittest.push_back(contestants[m]);
    }

    // 2.2) Erase the privileged contestants from the contestants list.
    contestants.erase(contestants.begin(), contestants.begin() + (pop_size / 10));

    // 3) Initiate the tournament.
    contestants_size = contestants.size();
    int i = 0, j = 0;

    do {

        do {

            i = rand() % contestants_size;
            j = rand() % contestants_size;

        // 3.1) FIXME: We're distinguishing contestants according to their 
        // evaluation function value. The impact of this may not be good.
        } while (contestants[i]->id == contestants[j]->id);

        if (contestants[i]->total_distance <= contestants[j]->total_distance) {

            fittest.push_back(contestants[i]);

            contestants.erase(
                    remove(
                            contestants.begin(), 
                            contestants.end(), 
                            contestants[i]), 
                    contestants.end());

        } else if (contestants[j]->total_distance < contestants[i]->total_distance) {

            fittest.push_back(contestants[j]);

            contestants.erase(
                    remove(
                            contestants.begin(), 
                            contestants.end(), 
                            contestants[j]), 
                    contestants.end());
        }

        contestants_size = contestants.size();

    } while (fittest.size() < (pop_size / 2));
}

//! Apply a mutation on a CVRP instance, i.e. a 1-interchange between routes 
//! and a 2-opt exchange in a random route.
/*!
*/
int Mutation(CVRP * graph) {

    // 1) Apply a 1-interchange procedure on the CVRP instance.
    graph->RandomOneInterchange();

    int stations_size = graph->stations.size();
    //int num_routes = graph->routes.size();
    int i = 0;

    // 2) Select 1 station at random from the CVRP's list, save a reference of 
    // its Route.
    i = rand() % (stations_size - 1) + 1;

    CVRP::Route * route_i = graph->stations[i]->route;
    graph->TwoOptExchange(route_i, CVRP::HARD_2_OPT);

    return 0;
}

//! Translates a TSP-like representation of a CVRP instance to an actual CVRP 
//! instance.
/*!
*/
int FromTSPRepresentation(
        CVRP * & offspring,
        CVRP::Station ** rep, 
        int req_num_stations,
        int req_num_routes,
        int rep_size) {

    //CVRP * translated = new CVRP();

    // 1) General graph attributes.
    long total_distance     = 0;

    // 1.1) FIXME: I find these things quite ugly...
    //offspring->capacity    = parent->capacity;

    // 1.2) Create a depot.
    CVRP::Station * depot = new CVRP::Station();

    // 1.3) Use the temporary offspring's depot to create a new dedicated depot.
    depot->index    = offspring->depot->index;
    depot->demand   = offspring->depot->demand;
    depot->coords   = offspring->depot->coords;
    depot->type     = offspring->depot->type;

    depot->next     = NULL;
    depot->prev     = NULL;
    depot->route    = NULL;

    offspring->depot = depot;

    // 1.3.1) Add the depot to the station's list.
    offspring->stations.push_back(depot);

    // 1.4) Remember that we should aim for a target number of routes, which 
    // may need to be re-adjusted given the TSP-like representation we've 
    // been given.
    unsigned int target_num_routes = req_num_routes;

    // 2) Cycle through the TSP-like representation array, and create Route and 
    // Station instances accordingly.
    int route_index         = 0;
    int route_qty_supplied  = 0;
    int route_num_stations  = 0;
    long route_distance     = 0;

    // 2.1) Pointers for working out new Stations.
    CVRP::Station * station         = NULL;
    CVRP::Station * prev_station    = NULL;

    // 2.2) Pointers for working out new Routes.
    CVRP::Route * curr_route        = NULL;

    // 2.3) Begin cycling the TSP-like representation.
    int i = 0;

    while (i < (rep_size - 1)) {

        // 2.3.1) Create a new Route.
        curr_route = new CVRP::Route();

        // 2.3.2) Leave the Routes alone for a little bit (actually just 
        // pick on them for just a little bit to update the start_station and 
        // end_station attributes) and handle the Stations now.
        prev_station = offspring->depot;

        while (rep[i]->type != CVRP::DEPOT) {

            station = new CVRP::Station();

            // 2.3.2.1) Copy the elements.
            station->index  = rep[i]->index;
            station->demand = rep[i]->demand;
            station->coords = rep[i]->coords;
            station->type   = rep[i]->type;
            station->route  = curr_route;

            // 2.3.2.2) Make the prev pointer of the new Station object to 
            // point to prev_station.
            station->prev   = prev_station;

            // 2.3.2.3) Make the next pointer of prev_station point to 
            // this new Station object.
            if (prev_station->type != CVRP::DEPOT) {

                prev_station->next = station;

            } else {

                curr_route->start_station = station;
            }

            // 2.3.2.4) If the next station of the graph being copied is the 
            // depot, make the new Station's next pointer point to its depot.
            if (rep[i + 1]->type == CVRP::DEPOT) {

                station->next = offspring->depot;

                curr_route->end_station = station;
            }

            // 2.3.2.5) Add the new object to the Station vector.
            offspring->stations.push_back(station);

            // 2.3.2.6) Update some variables for curr_route.
            route_qty_supplied += station->demand;
            route_num_stations++;
            route_distance += CVRP::Euc2DistCalc(station->prev, station);

            // 2.3.2.6) Advance.
            i++;
            prev_station = station;
        }

        // 2.3.3) Distance from end_station to the depot and update the 
        // graph's total distance.
        route_distance += CVRP::Euc2DistCalc(station, station->next);
        total_distance += route_distance;

        // 2.3.4) Update the remaining route's attributes.
        curr_route->index        = route_index++;
        curr_route->qty_supplied = route_qty_supplied;
        curr_route->num_stations = route_num_stations;
        curr_route->distance     = route_distance;
        curr_route->capacity     = offspring->capacity;
        curr_route->active       = true;

        // 2.3.5) Add the new route to the CVRP's Route list.
        offspring->routes.push_back(curr_route);

        // 2.3.6) Reset and advance.
        route_qty_supplied = 0;
        route_distance = 0;
        route_num_stations = 0;
        i++;
    }

    // 3) Check if a re-adjustment in the number of Routes of offspring is 
    // necessary.

    std::list <CVRP::Route *>::iterator refugees;
    std::list <CVRP::Route *>::iterator host;

    CVRP::Station * refugee = NULL;
    CVRP::Station * sheltered = NULL;

    long distance_change = 0;

    while (offspring->routes.size() > target_num_routes) {

        // Get a reference to the last element on the Route list (one assumes 
        // the list isn't empty).
        refugees = --offspring->routes.end();

        // The first refugee to allocate to a host is the first station on 
        // refugees.
        refugee = (*refugees)->start_station;

        for (
                host = offspring->routes.begin(); 
                host != refugees; 
                ) {

            // Go ahead with an insertion if the refugee Station does not 
            // exceed the capacity of the the host route 
            while (
                    refugee->type != CVRP::DEPOT 
                    && ((*host)->qty_supplied + refugee->demand < offspring->capacity)) {

                sheltered = refugee;
                refugee = refugee->next;

                (*host)->end_station->next = sheltered;

                sheltered->prev = (*host)->end_station;
                sheltered->next = offspring->depot;

                sheltered->route = (*host);

                // Update the info of the hosting route.
                (*host)->qty_supplied += sheltered->demand;

                distance_change = 
                        - CVRP::Euc2DistCalc((*host)->end_station, offspring->depot) +
                        + CVRP::Euc2DistCalc((*host)->end_station, sheltered) 
                        + CVRP::Euc2DistCalc(sheltered, sheltered->next);

                (*host)->distance += distance_change;
                total_distance += distance_change;

                // Update the host route's end station.
                (*host)->end_station = sheltered;

                (*host)->num_stations++;

                // Update the number of stations in the refugees Route.
                (*refugees)->num_stations--;
            }

            if (refugee->type == CVRP::DEPOT) {

                break;

            } else {

                host++;
            }
        }

        // If the refugees Route is empty, delete it, otherwise, abort the 
        // whole procedure.
        if ((*refugees)->num_stations > 0) {

            return -1;

        } else {

            // Update the total distance of the graph.
            total_distance = total_distance - (*refugees)->distance;

            delete(*refugees);
            offspring->routes.erase(refugees);
        }
    }

    // 3.x) 
    if ((int) offspring->routes.size() != req_num_routes
            || (int) offspring->stations.size() != req_num_stations) {

        return -1;
    }

    // 4) (Other) general graph attributes.
    offspring->num_stations    = offspring->stations.size();
    offspring->num_routes      = offspring->routes.size();
    offspring->total_distance  = total_distance;

    return 0;
}

//! Generate a TSP-like representation of a CVRP intance, for application of 
//! crossover operators such as PMX or OX.
/*!
*/
CVRP::Station ** CreateTSPRepresentation(CVRP * graph, int * & inter_depots) {

    int num_stations = graph->stations.size();
    int num_routes = graph->routes.size();

    // 1) Create an array of Station instances, with (num_routes) depots 
    // in between Routes, i.e. for a CVRP with 8 stations + depot D and 3 routes 
    // we want an array like this: [1][5][8][D][3][2][D][7][6][4][D].
    CVRP::Station ** path = (CVRP::Station **) calloc (
                                        num_stations + (num_routes), 
                                        sizeof(CVRP::Station*));

    // 2) Let us create the array in the form [1][5][8][D][3][2][D][7][6][4][D].
    std::list <CVRP::Route *>::iterator rt;
    inter_depots = (int *) calloc((num_routes), sizeof(int));
    int i = 0, j = 0;

    for (rt = graph->routes.begin(); rt != graph->routes.end(); rt++) {

        CVRP::Station * curr_station = (*rt)->start_station;

        // 2.1) Start copying the references, i.e. pointers, we're not 
        // creating new objects.
        while (curr_station->type != CVRP::DEPOT) {

            path[i++] = curr_station;
            curr_station = curr_station->next;
        }

        // 2.2) Leave some space for the intermediate depots.
        inter_depots[j++] = i++;
    }

    // 3) Set the intermediate Station between routes as the depot.
    for (int k = 0; k < (j - 1); k++) {

        path[inter_depots[k]] = graph->GetDepot();
    }

    // 4) Clean-up and return.
    free(inter_depots);
    return path;
}

//! Sort the route vector in a (capacity - qty_supplied) ascending 
//! fashion, i.e. lowest (capacity - qty_supplied) values first.
/*!
*/
bool SortRoutes(CVRP::Route * r_a, CVRP::Route * r_b) {

    return ((r_a->capacity - r_a->qty_supplied) 
            < (r_b->capacity - r_b->qty_supplied));
}

//! Generate 1 new CVRP instance via a crossover operation with genetic 
//! material from 2 parents. A variation of the Partially Mapped Crossover 
//! (PMX) method BRBAX (see Bermudez2010) is used.
/*!
*/
int CrossOver(
        CVRP * parent_i, 
        CVRP * parent_j, 
        CVRP * & cvrp_offspring, 
        double mutation_prob,
        int req_num_stations,
        int req_num_routes,
        int req_capacity) {

    // 1) Assuming that both parents have an equal number of stations and 
    // routes, any parent will do to keep track of the number of stations and 
    // routes.
    /*if (parent_i->capacity == 0 || parent_j->capacity == 0) {

        parent_i->PrintRoutes();
        parent_i->Validate();

        parent_j->PrintRoutes();
        parent_j->Validate();

        return -1;
    }*/

    int capacity        = req_capacity;
    int num_stations    = req_num_stations;
    int num_routes      = req_num_routes;

    // 1.1) The size of a TSP-like representation of a CVRP instance has this 
    // size, check CreateTSPRepresentation() to know why.
    //int tsp_rep_size    = num_stations + (num_routes);
    int tsp_rep_size    = 0;

    // 2) In order to apply PMX / BRBAX to the CVRP, the path representation is 
    // extended to contain multiple copies of the depot, each copy 
    // acting as a separator between routes (see section 6.5.3 on Toth2002).
    /*int * inter_depots_parent_1 = NULL;
    CVRP::Station ** parent_1 = CreateTSPRepresentation(
                                        parent_i,
                                        inter_depots_parent_1);

    int * inter_depots_parent_2 = NULL;
    CVRP::Station ** parent_2 = CreateTSPRepresentation(
                                        parent_j,
                                        inter_depots_parent_2);*/

    CVRP::Station ** offspring = (CVRP::Station **) calloc (
                                        num_stations + (2 * num_routes), 
                                        sizeof(CVRP::Station*));

    // 3) Start the crossover process.

    // 3.1) Find the best (num_routes / 2) routes of P1 with the lowest 
    // (capacity - qty_supplied) values.

    // FIXME: I use <std::list> for the graph->routes list due to iterator 
    // invalidation problems when adding / deleting CVRP::Route instances during 
    // the CWS algorithm. I now pass it to std::vector to use std::sort(). This 
    // is maybe inefficient, but I think it is necessary nonetheless.
    std::vector <CVRP::Route *> routes;

    // 3.1.1) Build up the route vector.
    std::list <CVRP::Route *>::iterator rt;

    for (rt = parent_i->routes.begin(); rt != parent_i->routes.end(); rt++) {

        routes.push_back(*rt);
    }

    // 3.1.2) Sort the route vector in a (capacity - qty_supplied) ascending 
    // fashion, i.e. lowest (capacity - qty_supplied) values first.
    sort(
            routes.begin(), 
            routes.end(), 
            SortRoutes);

    // 3.2) Array which marks the stations already existent on the offspring. 
    // If used_stations[i] == FALSE, the station with index i is available for 
    // insertion in the offspring.
    std::vector <bool> used_stations;

    // 3.2.1) As a starting point, all stations are available.
    for (int i = 0; i < num_stations; i++) {

        used_stations.push_back(false);
    }

    // 3.3) Pass the best (num_routes / 2) to the offspring.
    CVRP::Station * curr_station = NULL;

    // 3.3.1) For iterating over the TSP representation of the offspring. This 
    // will be used later on when filling up offspring[] from parent_2[].
    int j = 0;

    for (int i = 0; i < (num_routes / 2); i++) {

        curr_station = routes[i]->start_station;

        while (curr_station->type != CVRP::DEPOT) {

            offspring[j] = curr_station;

            // 3.3.2) Update used_stations.
            used_stations[offspring[j]->index] = true;

            j++;
            curr_station = curr_station->next;
        }

        // 3.3.3) Delimitation of a route in the TSP-like representation of a 
        // CVRP instance.
        offspring[j++] = parent_i->GetDepot();
    }

    // 3.4) Fill up the remaining positions of the offspring's TSP 
    // representation by attempting to keep the order of parent_2.

    // 3.4.1) Track the new routes capacity as you go.
    int max_new_qty_supplied = 0;
    int new_qty_supplied = 0;

    for (rt = parent_j->routes.begin(); rt != parent_j->routes.end(); rt++) {

        curr_station = (*rt)->start_station;

        new_qty_supplied = 0;

        while (curr_station->type != CVRP::DEPOT) {

            // 3.4.1) Only add the station to offspring[] if the station is 
            // not in offspring[] already.
            if (!used_stations[curr_station->index]) {

                offspring[j] = curr_station;

                // 3.4.2) Update used_stations. Not that it matters anyway...
                used_stations[offspring[j]->index] = true;

                // 3.4.3) Update new_qty_supplied.
                new_qty_supplied += offspring[j]->demand;

                j++;
            }

            curr_station = curr_station->next;
        }

        // 3.4.4) Update max_new_qty_supplied.
        if (max_new_qty_supplied < new_qty_supplied) {

            max_new_qty_supplied = new_qty_supplied;
        }

        // 3.4.5) Again, the route delimitation...
        if (offspring[j - 1]->type != CVRP::DEPOT)
            offspring[j++] = parent_j->GetDepot();
    }

    tsp_rep_size = j;

    // 3.5) The PMX / BRBAX crossover process is over. We might have ended up 
    // with a invalid CVRP instance (e.g. if route->qty_supplied exceeds the 
    // CVRP->capacity). Test this first before applying the TSP-to-CVRP 
    // translation.

    // 3.5.1) We wont be needed these TSP-like CVRP representations anymore.
    /*free(parent_1);
    free(parent_2);*/

    // 3.5.2) Just for fun, print the TSP-like CVRP representation.

    /*printf("OFFSPRING'S TSP REPRESENTATION:\n\t|");

    for (int m = 0; m < tsp_rep_size; m++) {

        printf("%3d|", offspring[m]->index - 1);
    }

    printf("\n\n");*/

    // 3.5.3) Abort (the word 'abort' has a not-so-nice double meaning in the 
    // crossover process used in genetic algorithms, right?) the offspring 
    // generation process (deleting all allocated memory) if the qty_supplied 
    // value of the new routes increases above capacity.
    if (max_new_qty_supplied > capacity) {

        /*printf("CROSSOVER: Offspring violates capacity constraint (%d > %d). "\
                "ABORTING.\n", 
                max_new_qty_supplied, 
                capacity);*/

        free(offspring);

        return -1;
    }

    // 3.6) It's OK to translate the TSP-like representation into a CVRP 
    // representation.
    cvrp_offspring->capacity = parent_i->capacity;
    cvrp_offspring->dataset = parent_i->dataset;

    // 3.6.1) Temporarily make the offspring's depot point to the parent's 
    // depot.
    cvrp_offspring->depot = parent_i->depot;

    // 3.6.2) Temporarily save the required number of routes in cvrp_offspring. 
    // Remember that we should aim for a target number of routes, which 
    // may need to be re-adjusted given the TSP-like representation we're 
    // passing to FromTSPRepresentation().
    cvrp_offspring->num_routes = parent_i->routes.size();

    if (FromTSPRepresentation(
            cvrp_offspring,
            offspring, 
            num_stations,
            num_routes,
            tsp_rep_size) < 0) {

        /*printf("CROSSOVER: Offspring violates route number constraint (%d > %d). "\
                    "ABORTING.\n", 
                (int) cvrp_offspring->routes.size(), 
                (int) parent_i->routes.size());*/

        free(offspring);

        return -1;
    }

    free(offspring);

    // 3.7) According to the probability of mutations, apply one to the 
    // offspring (this only affects route distance, so no infeasible solutions 
    // are generated here).
    double r = (double) rand() / (double) (RAND_MAX + 1.0);

    if (mutation_prob > r) {

        // 3.7.1) (Un)lucky offspring. Considering the type of mutation, I 
        // would say it's most probably unlucky (route(s) will get crossings).
        Mutation(cvrp_offspring);
    }

    //cvrp_offspring->Validate();

    cvrp_offspring->id = global_id++;

    return 0;
}

int GenerateRandomIndividual(
        CVRP * & individual,
        std::string dataset, 
        int num_routes) {

    individual = new CVRP(dataset);

    std::vector <int> station_index_pool;
    std::vector <int> route_index_pool;

    int station_index_pool_size = individual->stations.size() - 1;

    for (int i = 0; i < station_index_pool_size; i++) {

        station_index_pool.push_back(i + 1);
    }

    int route_index_pool_size = num_routes;

    for (int i = 0; i < route_index_pool_size; i++) {

        route_index_pool.push_back(i);
    }

    // Pointers for current stations and routes when handled within cycles.
    CVRP::Station * station = NULL;
    CVRP::Route * route     = NULL;

    // Even though the CVRP instance saves routes in a std::list, let's use 
    // a vector for the building phase, since we'll need random access to each
    // route.
    std::vector <CVRP::Route *> routes;

    for (int i = 0; i < num_routes; i++) {

        // Create a new route to get stations into.
        route = new CVRP::Route();

        // Handle the initial attributes of the new route.
        route->index        = i;
        route->qty_supplied = 0;
        route->capacity     = individual->capacity;
        route->num_stations = 0;
        route->distance     = 0;
        route->active       = true;

        // Initially set the start and end pointers of the routes to the depot.
        route->start_station    = individual->depot;
        route->end_station      = individual->depot;

        // Add the route to the route vector.
        routes.push_back(route);
    }

    long distance_change    = 0;
    long total_distance     = 0;

    int i = 0, j = 0, station_index = 0, route_index = 0;

    // FIXME: I'm trusting that an eternal cycle in which a station can never 
    // find a route to get in will never happen... Maybe I should not be so 
    // optimistic.
    while (station_index_pool_size > 0) {

        i = rand() % station_index_pool_size;
        station_index = station_index_pool[i];

        station = individual->stations[station_index];

        /*printf("GENERATE RANDOM: Picked %d, station[%d] with index %d.\n",
                i,
                station_index_pool[i],
                individual->stations[station_index]->index - 1);*/

        // Get a random route index.
        j = rand() % route_index_pool_size;
        route_index = route_index_pool[j];

        /*printf("GENERATE RANDOM: Picked %d, route[%d] with index %d.\n",
                j,
                route_index_pool[j],
                routes[route_index]->index);*/

        // Keep looking for a suitable route index until the vehicle capacity 
        // constraint is not violated.
        while (routes[route_index]->qty_supplied + station->demand > individual->capacity) {

            /*printf("GENERATE RANDOM: Route %d, route[%d] with index %d is full.\n",
                    j,
                    route_index_pool[j],
                    routes[route_index]->index);*/

            // Erase the unfeasible index from the route index pool, so that 
            // one doesn't choose that again.
            //route_index_pool.erase(route_index_pool.begin() + j);
            route_index_pool.erase(
                    remove(
                            route_index_pool.begin(), 
                            route_index_pool.end(), 
                            route_index_pool[j]), 
                    route_index_pool.end());
            route_index_pool_size = route_index_pool.size();

            if (!(route_index_pool_size > 0)) {

                printf("GENERATE RANDOM: Individual violates route number "\
                            "constraint. ABORTING.\n");

                return -1;
            }

            j = rand() % route_index_pool_size;
            route_index = route_index_pool[j];

            /*printf("GENERATE RANDOM: Picked %d, route[%d] with index %d.\n",
                    j,
                    route_index_pool[j],
                    routes[route_index]->index);*/
        }

        // Add the station to the route, and handle all necessary attributes of 
        // Station and Route.
        station->route = routes[route_index];

        // If this is the first station to be added to a route, set it as the 
        // start station.
        if (routes[route_index]->start_station->type == CVRP::DEPOT) {

            routes[route_index]->start_station = station;
        }

        station->next = individual->depot;
        station->prev = routes[route_index]->end_station;

        // Also, update the previous stations' next pointer.
        if (station->prev->type != CVRP::DEPOT) {

            station->prev->next = station;
        }

        // Update the quantitative route's attributes.
        distance_change = 
                - CVRP::Euc2DistCalc(routes[route_index]->end_station, individual->depot) +
                + CVRP::Euc2DistCalc(routes[route_index]->end_station, station) 
                + CVRP::Euc2DistCalc(station, station->next);

        routes[route_index]->distance += distance_change;
        routes[route_index]->qty_supplied += station->demand;
        routes[route_index]->num_stations++;

        // Finally update the route's end_station pointer.
        routes[route_index]->end_station = station;

        // Erase the index of the freshly added Station from the station index 
        // pool, so that one doesn't choose that again.
        //station_index_pool.erase(station_index_pool.begin() + i);
        station_index_pool.erase(
                remove(
                        station_index_pool.begin(), 
                        station_index_pool.end(), 
                        station_index_pool[i]), 
                station_index_pool.end());

        station_index_pool_size = station_index_pool.size();

        /*printf("NEW INDIVIDUAL'S INDEX POOL:\n\t|");

        for (int s = 0; s < station_index_pool_size; s++) {

            printf("[%3d]->[%3d]|", station_index_pool[s], individual->stations[station_index_pool[s]]->index - 1);
        }

        printf("\n");*/
    }

    // All stations have been attributed a route. Let's wrap everything up...

    // Create the individual's CVRP route list.
    for (int k = 0; k < num_routes; k++) {

        individual->routes.push_back(routes[k]);
        total_distance += routes[k]->distance;
    }

    // Update the few remaining attributes of individual.
    individual->num_stations    = individual->stations.size();
    individual->num_routes      = individual->routes.size();
    individual->total_distance  = total_distance;

    return 0;
}

void GenerateInitialPopulation(
        CVRP * cws_parallel,
        CVRP * cws_sequential,
        CVRP * sa,
        long pop_size,
        std::vector <CVRP *> & population) {

    // 1) Add the CWS and SA solutions to the population.
    //population.push_back(cws_parallel);
    //population.push_back(cws_sequential);
    //population.push_back(sa);

    // 2) Generate the remaining population elements.
    int i = 0;

    CVRP * new_individual = NULL;

    while (i < (pop_size)) {

        if (!(GenerateRandomIndividual(
                new_individual,
                cws_parallel->dataset,
                (int) cws_parallel->routes.size()) < 0)) {

            new_individual->id = global_id++;
            //new_individual->PrintRoutes();
            //new_individual->Validate();
            population.push_back(new_individual);

            i++;
        }
    }
}

void Evolutionary::Run(
        CVRP * cws_parallel,
        CVRP * cws_sequential,
        CVRP * sa,
        uint32_t max_generations,
        long pop_size,
        double mutation_prob) {

    uint32_t generation = 0;
    uint32_t S_b_generation = 0;

    // 1) Generate the initial population, including the solutions generated by
    // the CWS algorithm, and those obtained via SA from multiple parameters.
    // FIXME: In the case of the CWS algorithm solutions, both must have the 
    // same number of routes for this to work properly.
    std::vector <CVRP *> population;
    std::vector <CVRP *> fittest;
    std::vector <CVRP *> new_generation;

    // 1.1) FIXME: This should be done elsewhere, I know, but let's go with this 
    // for now.
    cws_parallel->id = global_id++;
    cws_sequential->id = global_id++;
    sa->id = global_id++;

    GenerateInitialPopulation(
            cws_parallel,
            cws_sequential,
            sa,
            pop_size,
            population);

    // 1.2) Save the best solution S_b.
    sort(
        population.begin(), 
        population.end(), 
        SortCandidates);

    CVRP * S_b = new CVRP();
    S_b->CopyFrom(population.front());
    long initial_best = S_b->total_distance;

    // 1.3) Start counting the execution time of the algorithm
    struct timeval tm;
    gettimeofday(&tm, NULL);
    double start_time = (tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0);

    // 1.4) Saving some hard constraints for the CrossOver() procedure, so that 
    // appropriate CVRP instances are generated.
    int req_num_stations = cws_parallel->stations.size();
    int req_num_routes = cws_parallel->routes.size();
    int req_capacity = cws_parallel->capacity;

    // 2) Initialize the evolutionary cycle.
    do {

        // 2.1) Random tournament selection, until the number of selected 
        // winners (or contestants) is pop_size / 2. Notice that parents and 
        // their offspring fight between each other, and both are possibly 
        // chosen.
        RandomTournamentSelection(
                population,
                fittest);

        // 2.2) Take care of the 'dead'.
        std::vector <CVRP *>::iterator it;

        for (it = population.begin(); it != population.end(); it++) {
            delete (*it);
        }

        population.clear();

        // 2.3) Apply CROSSOVER and MUTATION to generate offspring. Note that 
        // offspring is created from pairs of solutions. So, (pop_size / 2) 
        // pairs are selected at random, in order to generate (pop_size / 2) 
        // offspring, so that pop_size remains constant.
        int i = 0, j = 0;
        int fittest_size = fittest.size();

        CVRP * offspring = NULL;

        while (new_generation.size() < (pop_size / 2)) {

            do {

                i = rand() % fittest_size;
                j = rand() % fittest_size;

            } while (fittest[i]->id == fittest[j]->id);

            /*printf("EVOLUTIONARY: Chose 2 fittest %d and %d\n",
                    i, j);*/

            /*fittest[i]->PrintRoutes();
            fittest[i]->Validate();

            fittest[j]->PrintRoutes();
            fittest[j]->Validate();*/

            offspring = new CVRP();

            if (!(CrossOver(
                    fittest[i], 
                    fittest[j], 
                    offspring,
                    mutation_prob,
                    req_num_stations,
                    req_num_routes, 
                    req_capacity) < 0)) {

                new_generation.push_back(offspring);

            } else {

                delete offspring;
            }
        }

        // 2.4) Join the parents and offspring into a unified population.
        int new_generation_size = new_generation.size();

        for (int i = 0; i < fittest_size; i++) {

            population.push_back(fittest[i]);
        }

        fittest.clear();

        for (int i = 0; i < new_generation_size; i++) {

            population.push_back(new_generation[i]);
        }

        new_generation.clear();

        // 2.5) Save the best solution S_b so far.
        sort(
                population.begin(), 
                population.end(), 
                SortCandidates);

        if (S_b->total_distance > population.front()->total_distance) {

            S_b->CopyFrom(population.front());
            S_b_generation = generation;
        }

        printf("EVOLUTIONARY: At generation %d: "\
                "\n\t[S_b] = %ld %d"\
                "\n\t[S_c] = %ld %d"\
                "\n\t[P] = %d\n\n",
                generation,
                S_b->total_distance, (int) S_b->stations.size(),
                population.front()->total_distance, (int) population.front()->stations.size(),
                (int) population.size());

    } while (generation++ < max_generations);

    gettimeofday(&tm, NULL);
    double ga_time = ((tm.tv_sec * 1000.0) + (tm.tv_usec/1000.0)) - start_time;

    // Apply a 2-opt de-crossing procedure over the several routes.
    std::list <CVRP::Route *>::iterator rt;

    for (rt = S_b->routes.begin(); rt != S_b->routes.end(); rt++) {

        S_b->TwoOptExchange((*rt), CVRP::HARD_2_OPT);
    }

    // Print a Evolutionary report.
    PrintEvolutionaryReport(
            S_b->dataset.c_str(), 
            initial_best,
            S_b->total_distance,
            S_b_generation,
            pop_size,
            mutation_prob, 
            generation, 
            ga_time);

    S_b->PrintRoutes();
    S_b->Validate();
}

