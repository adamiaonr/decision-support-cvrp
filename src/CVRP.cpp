#include "CVRP.h"
#include "DatasetParser.h"

//! CVRP constructor.
/*!
*/
CVRP::CVRP()
{
    this->dataset           = std::string("UNKNOWN");
    this->num_stations      = 0;
    this->num_routes        = 0;
    this->total_distance    = 0;
    this->capacity          = 0;

    this->depot = NULL;
}

//! Alternative CVRP constructor, sets the dataset file name and fills the CVRP 
//! instance from a dataset file in the format as in format as in 
//! http://www.coin-or.org/SYMPHONY/branchandcut/VRP/index.htm.
/*!
*/
CVRP::CVRP(std::string dataset)
{
    this->dataset           = dataset;
    this->num_stations      = 0;
    this->num_routes        = 0;
    this->total_distance    = 0;
    this->capacity          = 0;

    this->depot = NULL;

    if (this->FillFromDataset(dataset.c_str()) < 0) {

        printf("CVRP CREATION: Error while filling CVRP instance from dataset "\
                "%s.\n", 
                dataset.c_str());
    }
}

//! CVRP destructor.
/*!
*/
CVRP::~CVRP()
{
    this->num_stations      = 0;
    this->num_routes        = 0;
    this->total_distance    = 0;
    this->capacity          = 0;

    this->depot = NULL;

    // 1) Erase the previous Station elements.
    std::vector <Station *>::iterator it;

    for (it = this->stations.begin(); it != this->stations.end(); it++) {

        delete (*it);
    }

    this->stations.clear();

    // 2) Erase the previous Route elements.
    std::list <Route *>::iterator rt;

    for (rt = this->routes.begin(); rt != this->routes.end(); ) {

        delete (*rt);
        rt = this->routes.erase(rt);
    }
}

//! Fills the CVRP instance from a dataset file in the format as in format as 
//! in http://www.coin-or.org/SYMPHONY/branchandcut/VRP/index.htm.
/*!
*/
int CVRP::FillFromDataset(const char * dataset) {

    return DatasetParser::ParseDataset(dataset, this);
}

//! Copies the contents of the CVRP * given as argument into the calling 
//! CVRP object.
/*!
*/
void CVRP::CopyFrom(CVRP * graph) {

    // General graph attributes.
    this->dataset           = graph->dataset;
    this->num_stations      = graph->num_stations;
    this->num_routes        = graph->num_routes;
    this->total_distance    = graph->total_distance;
    this->capacity          = graph->capacity;

    // Create a Station object representing the depot.
    Station * depot = new Station();

    depot->index    = graph->depot->index;
    depot->demand   = graph->depot->demand;
    depot->coords   = graph->depot->coords;
    depot->type     = graph->depot->type;

    // All the remaining elements are NULL
    depot->next     = NULL;
    depot->prev     = NULL;
    depot->route    = NULL;

    this->depot = depot;

    // For the lists, let's actually create copies of the Station and Route 
    // objects pointed by list elements.

    // 1) Erase the previous Station elements.
    std::vector <Station *>::iterator it;

    for (it = this->stations.begin(); it != this->stations.end(); it++) {

        delete (*it);
    }

    this->stations.clear();

    // 2) Erase the previous Route elements.
    std::list <Route *>::iterator rt;

    for (rt = this->routes.begin(); rt != this->routes.end(); ) {

        delete (*rt);
        rt = this->routes.erase(rt);
    }

    this->stations.push_back(this->depot);

    // 3) Start by copying the Routes and copy the Station objects as you 
    // cycle through each route.
    for (rt = graph->routes.begin(); rt != graph->routes.end(); rt++) {

        Route * route = new Route();

        route->index        = (*rt)->index;
        route->qty_supplied = (*rt)->qty_supplied;
        route->capacity     = (*rt)->capacity;
        route->num_stations = (*rt)->num_stations;
        route->distance     = (*rt)->distance;
        route->active       = (*rt)->active;

        Station * curr_station = (*rt)->start_station;
        Station * prev_station = this->depot;
        Station * station;

        while (curr_station->type != CVRP::DEPOT) {

            station = new Station();

            // Copy the elements.
            station->index  = (curr_station)->index;
            station->demand = (curr_station)->demand;
            station->coords = (curr_station)->coords;
            station->type   = (curr_station)->type;
            station->route  = route;

            // Make the prev pointer of the new Station object to 
            // point to prev_station.
            station->prev   = prev_station;

            // Make the next pointer of prev_station point to 
            // this new Station object.
            if (prev_station->type != CVRP::DEPOT) {

                prev_station->next = station;

            } else {

                route->start_station = station;
            }

            // If the next station of the graph being copied is the 
            // depot, make the new Station's next pointer point to its depot.
            if (curr_station->next->type == CVRP::DEPOT) {

                station->next = this->depot;

                route->end_station = station;
            }

            // Add the new object to the Station vector.
            this->stations.push_back(station);

            // Advance.
            prev_station = station;
            curr_station = curr_station->next;
        }

        // Add the new route to the Route list.
        this->routes.push_back(route);
    }
}

long CVRP::Euc2DistCalc(euc_2d_coords i, euc_2d_coords j) {

    long xd = i.x - j.x;
    long yd = i.y - j.y;

    return (long) ceil(sqrt((xd * xd) + (yd * yd)));
}

long CVRP::Euc2DistCalc(Station * i, Station * j) {

    long xd = i->coords.x - j->coords.x;
    long yd = i->coords.y - j->coords.y;

    return (long) ceil(sqrt((xd * xd) + (yd * yd)));
}

//! Adds a Station (vertex) to a CVRP's graph. Sets the Station 
//! properties according to the parameters passed to the function.
/*!
*/
CVRP::Station * CVRP::AddStation(
        int index,
        int demand,
        euc_2d_coords coords,
        bool is_depot) {

    Station * new_station = new Station;

    if (new_station == NULL) {

        return NULL;
    }

    new_station->index = index;
    new_station->demand = demand;
    new_station->coords = coords;

    if (is_depot)
        new_station->type = CVRP::DEPOT;
    else
        new_station->type = CVRP::NORMAL;

    // Add a new entry to the CVRP's graph station vector.
    this->stations.push_back(new_station);

    // Update the CVRP's graph properties.
    this->num_stations++;

    if (new_station->type == CVRP::DEPOT)
        this->depot = new_station;

    return new_station;
}

CVRP::Station * CVRP::GetNextStation(Station * station) {

    return station->next;
}

CVRP::Station * CVRP::GetPrevStation(Station * station) {

    return station->prev;
}

void CVRP::LinkStations(Station * start, Station * end) {

    // If the start/end Station isn't the depot, check if there are adjacent 
    // vertices first. If that's the case, remove them, then proceed with the 
    // new edge addition (i.e. linking operation).

    // FIXME: Although this now seems kind of unproblematic now, one should 
    // check if this causes problems in the future.
    if (start->type != DEPOT && end->type != DEPOT) {

        Station * next_station = start->next; next_station->prev = NULL;
        Station * previous_station = end->prev; previous_station->next = NULL;
    }

    // Establish the link... That's it!
    start->next = end;
    end->prev = start;
}

/*void CVRP::DeleteStation(Station * station) {

    // Remove the Station * from the stations vector.
    this->stations.remove(station);

    // Delete the station instance;
    delete station;
}*/

void CVRP::PrintStations() {

    std::vector <Station *>::iterator it;

    printf("STATION LIST:\n");

    for (it = this->stations.begin(); it != this->stations.end(); it++) {

        printf("\t[INDEX]: %3d, [DEMAND]: %3d [COORDS]: {%3d,%3d}", 
                (*it)->index, 
                (*it)->demand, 
                (*it)->coords.x, 
                (*it)->coords.y);

        if ((*it)->type == DEPOT) {

            printf(", [DEPOT]: Y");
        }

        printf("\n");
    }
}

CVRP::Station * CVRP::GetDepot() {

    return this->depot;
}

CVRP::Route * CVRP::AddRoute(
        int qty_supplied,
        long distance,
        Station * start_station,
        Station * end_station) {

    Route * new_route = new Route;

    new_route->index = this->routes.size();
    new_route->qty_supplied = qty_supplied;
    new_route->capacity = this->capacity;
    new_route->distance = distance;
    new_route->num_stations = 0;

    // Also, update the total graph's distance.
    this->total_distance += distance;

    new_route->start_station  = start_station;
    new_route->end_station  = end_station;

    this->routes.push_back(new_route);

    Station * i = start_station;

    while (i->type != CVRP::DEPOT) {

        new_route->num_stations++;
        i = i->next;
    }

    return new_route;
}

void CVRP::DeleteRoute(Route * route) {

    // Remove the Route * from the CVRP's graph routes list.
    this->routes.remove(route);

    // Delete the route object;
    delete route;
}

CVRP::Route * CVRP::Assign2Route(Station * station, Route * route) {

    station->route = route;

    return station->route;
}

void CVRP::PrintRoutes() {

    printf("ROUTE LIST:\n");

    int k = 0;
    std::list <Route *>::iterator it;

    // Go through the list of current routes held by this graph. 
    // For each route, go from a Station u to the adjacent Station v (only 1 in 
    // this graph, considering that u != depot), printing its index as you go.
    for (it = this->routes.begin(); it != this->routes.end(); it++) {

        // Vertex descriptor for the initial station (i.e. immediately after 
        // the depot) of the route.
        Station * current_station = (*it)->start_station;

        printf("\n\tROUTE #%2d: %3d", 
                //(*it)->GetIndex(),
                ++k,
                (this->depot->index - 1));

        do {

            printf("|%3d", current_station->index - 1);
            current_station = current_station->next;

        } while (current_station->type != CVRP::DEPOT);

        printf("|%3d,  [D:%3ld][QS:%3d][NS:%3d]\n",
                (this->depot->index - 1),
                (*it)->distance,
                (*it)->qty_supplied,
                (*it)->num_stations);
    }
}

long CVRP::SwapStations(Station * i, Station * j, int num_stations) {

    Station * prev_i, * next_j, * aux, * pos_i, * pos_j, * k, * l;

    // Assuming that Stations i and j belong to the same route, recalculate 
    // the route's distance and return it.

    long new_distance = 
            + i->route->distance 
            + Euc2DistCalc(i->prev, j)
            + Euc2DistCalc(i, j->next)
            - Euc2DistCalc(i->prev, i)
            - Euc2DistCalc(j, j->next);

    // First, update the stations which are adjacent to the left and right of 
    // i and j.
    prev_i = i->prev;
    next_j = j->next;

    prev_i->next = j;
    next_j->prev = i;

    // Update the route's start or end station, if necessary.
    if (prev_i->type == CVRP::DEPOT) {

        i->route->start_station = j;

    } else if (next_j->type == CVRP::DEPOT) {

        i->route->end_station = i;
    }

    /*printf("SWAP: Left bound %3d now pointing to %3d (%3d)\n", 
            i->prev->index, i->prev->next->index, j->index);

    printf("SWAP: Right bound %3d now pointing to %3d (%3d)\n", 
            j->next->index, j->next->prev->index, i->index);*/

    // Second, the stations to be swapped.
    pos_i = i->next;
    pos_j = j->prev;

    i->prev = pos_i;
    i->next = next_j;

    j->next = pos_j;
    j->prev = prev_i;

    /*printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
            i->index, i->next->index, i->index, i->prev->index);

    printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
            j->index, j->next->index, j->index, j->prev->index);*/

    int m = 0;

    while (m++ < ((num_stations / 2))) {

        // Save the Station * to where one should go next.
        k = pos_i->next;
        l = pos_j->prev;

        // Swap pointers which point to the inside of the route.
        aux = pos_i->prev;
        pos_i->prev = pos_i->next;
        pos_i->next = aux;

        /*printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
                pos_i->index, pos_i->next->index, pos_i->index, pos_i->prev->index);*/

        if (pos_i == pos_j)
            break;

        aux = pos_j->next;
        pos_j->next = pos_j->prev;
        pos_j->prev = aux;

        /*printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
                pos_j->index, pos_j->next->index, pos_j->index, pos_j->prev->index);*/

        // Advance along the route.
        pos_i = k;
        pos_j = l;
    }

    // Update the graph's and route's distances.
    this->total_distance += (new_distance - i->route->distance);
    i->route->distance = new_distance;

    return i->route->distance;
}

CVRP::Station** CVRP::CreateRouteRep(CVRP::Route * route) {

    CVRP::Station** route_stations = (CVRP::Station**) calloc (
                                                    route->num_stations + 2, 
                                                    sizeof(CVRP::Station*));

    CVRP::Station * curr_route_station = route->start_station;

    // Just a reference to the graph's depot.
    CVRP::Station * depot = this->GetDepot();

    for (int i = 1; i < route->num_stations + 2; i++) {

        route_stations[i] = curr_route_station;
        curr_route_station = curr_route_station->next;
    }

    // Set the edges of the array to point to the depot.
    route_stations[0] = depot;
    route_stations[route->num_stations + 1] = depot;

    return route_stations;
}

//! Inserts a Station i (from Route a route_i) in a Route route_j, after 
//! Station j.
/*!
*/
long CVRP::InsertIntoRoute(
        CVRP::Station * i,
        CVRP::Route * route_j,
        CVRP::Station * j) {

    // P.1) If i is the depot, abort.
    if (i->type == CVRP::DEPOT) {

        printf("INSERT: Station[%d] is the depot. ABORTING.\n", i->index);

        return -1;
    }

    // P.2) If j is the depot, abort.
    if (j->type == CVRP::DEPOT) {

        printf("INSERT: Station[%d] is the depot. ABORTING.\n", j->index);

        return -1;
    }

    // P.3) Also, let's avoid depleting route_i.
    Route * route_i = i->route;

    if (route_i->num_stations < 2) {

        printf("INSERT: Route[%d] is too small. ABORTING.\n", i->route->index);

        return -1;
    }

    // 1.1) Pointers to the Stations [i - 1] and [i + 1].
    Station * before_i = i->prev;
    Station * after_i = i->next;

    // 1.2) Pointer to the Station [j + 1].
    Station * after_j = j->next;

    // 2.1) Update the pointers of elements [j] and [j + 1], in route_j.
    j->next = i;

    if (after_j->type != CVRP::DEPOT)
        after_j->prev = i;

    // 2.2) Update the pointers of elements [i - 1] and [i + 1], in route_i.
    if (before_i->type != CVRP::DEPOT)
        before_i->next = after_i;

    if (after_i->type != CVRP::DEPOT)
        after_i->prev = before_i;

    // 2.3) Update the prev and next pointers of i.
    i->prev = j;
    i->next = after_j;

    // 3.1) Update the route of i.
    i->route = route_j;

    // 3.2) If i is the new starting point of route_j, update it.
    if (i->prev->type == CVRP::DEPOT)
        i->route->start_station = i;

    // 3.3) Similarly, if the last station...
    if (i->next->type == CVRP::DEPOT)
        i->route->end_station = i;

    // 3.4) Do the same with route_i.
    if (before_i->type == CVRP::DEPOT)
        route_i->start_station = after_i;

    if (after_i->type == CVRP::DEPOT)
        route_i->end_station = before_i;

    // 4) Update distances.
    /*long new_j_distance =
            + route_j->distance
            + CVRP::Euc2DistCalc(j, i)
            + CVRP::Euc2DistCalc(i, after_j)
            - CVRP::Euc2DistCalc(j, after_j);

    long new_i_distance =
            + route_i->distance
            + CVRP::Euc2DistCalc(before_i, after_i)
            - CVRP::Euc2DistCalc(before_i, i)
            - CVRP::Euc2DistCalc(i, after_i);*/
    long new_j_distance = this->CalcRouteDist(route_j);
    long new_i_distance = this->CalcRouteDist(route_i);

    /*printf("INSERT: i(%ld,%ld) j(%ld,%ld) total(%ld)\n", 
            route_i->distance, new_i_distance,
            route_j->distance, new_j_distance,
            this->total_distance);*/

    this->total_distance += 
        + new_i_distance
        + new_j_distance
        - route_i->distance
        - route_j->distance;

    /*printf("INSERT: TOTAL DISTANCE BEFORE 2-OPT: %ld\n", 
            this->total_distance);*/

    // 4.1) Route j
    route_j->distance = new_j_distance;

    // 4.2) Route i
    route_i->distance = new_i_distance;

    // 4.4) Also, update the demand quantities served in each route.
    i->route->qty_supplied += i->demand;
    route_i->qty_supplied -= i->demand;

    // 4.5) ... aaaaand the number of stations.
    i->route->num_stations++;
    route_i->num_stations--;

    //this->PrintRoutes();

    return this->total_distance;
}

long CVRP::SwapStations(
        CVRP::Station * i,
        CVRP::Route * route_i,
        CVRP::Station * j,
        CVRP::Route * route_j) {

    Station * before_i = i->prev;
    Station * after_i = i->next;

    Station * before_j = j->prev;
    Station * after_j = j->next;

    /*printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
            i->index, i->next->index, i->index, i->prev->index);

    printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
            j->index, j->next->index, j->index, j->prev->index);*/

    // 1) Update the pointers of the elements before and after i and j
    before_i->next = j;
    after_i->prev = j;

    before_j->next = i;
    after_j->prev = i;

    // 2) Update the prev and next pointers of i and j.
    i->prev = before_j;
    i->next = after_j;

    // 2.1) Update i's route to route_j.
    i->route = route_j;

    // 2.2) If i is the new starting point of route j, update it.
    if (i->prev->type == CVRP::DEPOT)
        i->route->start_station = i;

    // 2.3) Similarly, if the last station...
    if (i->next->type == CVRP::DEPOT)
        i->route->end_station = i;

    j->prev = before_i;
    j->next = after_i;

    j->route = route_i;

    if (j->prev->type == CVRP::DEPOT)
        j->route->start_station = j;

    if (j->next->type == CVRP::DEPOT)
        j->route->end_station = j;

    /*printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
            i->index, i->next->index, i->index, i->prev->index);

    printf("SWAP: %3d->next = %3d %3d->prev = %3d\n", 
            j->index, j->next->index, j->index, j->prev->index);*/

    // 3) New distances for route_i and route_j.
    /*long new_i_distance =
        + route_i->distance
        + CVRP::Euc2DistCalc(before_i, j)
        + CVRP::Euc2DistCalc(j, after_i)
        - CVRP::Euc2DistCalc(before_j, j)
        - CVRP::Euc2DistCalc(j, after_j);

    long new_j_distance =
        + route_j->distance
        + CVRP::Euc2DistCalc(before_j, i)
        + CVRP::Euc2DistCalc(i, after_j)
        - CVRP::Euc2DistCalc(before_i, i)
        - CVRP::Euc2DistCalc(i, after_i);*/

    long new_j_distance = this->CalcRouteDist(route_j);
    long new_i_distance = this->CalcRouteDist(route_i);

    this->total_distance += 
        + new_i_distance
        + new_j_distance
        - route_i->distance
        - route_j->distance;

    /*printf("SWAP: TOTAL DISTANCE BEFORE 2-OPT: %ld\n", 
            this->total_distance);*/

    /*printf("SWAP: TOTAL DISTANCE: %ld \n\n", 
            this->total_distance);*/

    // 3.1) Route i
    route_i->distance = new_i_distance;

    // 3.2) Route j
    route_j->distance = new_j_distance;

    // 3.3) Also, update the demand quantities served in each route.
    i->route->qty_supplied += (i->demand - j->demand);
    j->route->qty_supplied += (j->demand - i->demand);

    return this->total_distance;
}

bool CVRP::Validate() {

    std::vector<bool> validation_array;
    bool validation = false;
    long distance = 0;

    for (int i = 0; i < this->num_stations; i++) {

        validation_array.push_back(true);
    }

    std::list <Route *>::iterator rt;

    for (rt = this->routes.begin(); rt != this->routes.end(); rt++) {

        Station * curr_station = (*rt)->start_station;
        Station * prev_station = NULL;

        while (curr_station->type != CVRP::DEPOT) {

            // Check if this position was already taken.
            if (validation_array[curr_station->index - 1] == false) {

                validation = true;
            }

            validation_array[curr_station->index - 1] = false;
            distance += CVRP::Euc2DistCalc(curr_station->prev, curr_station);

            // Advance.
            prev_station = curr_station;
            curr_station = curr_station->next;
        }

        distance += CVRP::Euc2DistCalc(prev_station, prev_station->next);

        validation_array[curr_station->index - 1] = false;
    }

    int validation_array_size = validation_array.size();

    for (int i = 0; i < validation_array_size; i++) {

        validation = (validation | validation_array[i]);

        /*printf("ROUTE VALIDATION: station[%d] %s in solution\n", 
                i,
                (validation_array[i] ? "NOT" : "IS"));*/
    }

    printf("\nSOLUTION VALIDATION REPORT: "\
            "\n\t[NUM STATIONS] = %d %d"\
            "\n\t[STATION INCLUSION] = %s\n\t[MEAS. TOTAL DISTANCE] = %ld (SAVED IS %ld)\n\n", 
            (int) this->num_stations, (int) this->stations.size(),
            (validation ? "NOT OK" : "OK"),
            distance,
            this->total_distance);

    return !(validation);
}

long CVRP::CalcRouteDist(CVRP::Route * rt) {

    long distance = 0;

    Station * curr_station = (rt)->start_station;

    while (curr_station->next->type != CVRP::DEPOT) {

        distance += CVRP::Euc2DistCalc(curr_station->prev, curr_station);

        // Advance.
        curr_station = curr_station->next;
    }

    distance += CVRP::Euc2DistCalc(curr_station->prev, curr_station);
    distance += CVRP::Euc2DistCalc(curr_station, curr_station->next);

    return distance;
}

long CVRP::OneInterchange(CVRP::IChange & ichange) {

    /*printf("1-INTERCHANGE W/ PARAMS: "\
            "\n\tstation[%d] = %d"\
            "\n\tstation[%d] = %d"\
            "\n\troute[%d] = %d"\
            "\n\troute[%d] = %d"\
            "\n\tOP = %d"\
            "\n\t[DISTANCE (OR.)] = %ld"\
            "\n\t[DISTANCE (NEW)] = %ld\n\n", 
            ichange.i, stations[ichange.i]->index - 1,
            ichange.j, stations[ichange.j]->index - 1,
            ichange.i, stations[ichange.i]->route->index,
            ichange.j, stations[ichange.j]->route->index,
            ichange.operation,
            this->total_distance,
            ichange.distance);*/

    Station * station_i = this->stations[ichange.i];
    Station * station_j = this->stations[ichange.j];

    Route * route_i = this->stations[ichange.i]->route;
    Route * route_j = this->stations[ichange.j]->route;

    int qty_supplied_i = 0;
    int qty_supplied_j = 0;

    int demand_i = 0;
    int demand_j = 0;

    long distance = 0;

    // Operator (0,1), i.e. element j goes into route of element i.
    if (ichange.operation == ICHANGE_OP01) {

        qty_supplied_i = station_i->route->qty_supplied;
        demand_j = station_j->demand;

        // First verify if the change is feasible.
        if (qty_supplied_i + demand_j <= this->capacity) {

            this->InsertIntoRoute(
                    station_j,
                    station_i->route,
                    station_i);

            this->TwoOptExchange(station_i->route, CVRP::HARD_2_OPT);
            this->TwoOptExchange(station_j->route, CVRP::HARD_2_OPT);

            distance = this->total_distance;
        }
    // Operator (1,0), i.e. element i goes into route of element j.
    } else if (ichange.operation == ICHANGE_OP10) {

        qty_supplied_j = station_j->route->qty_supplied;
        demand_i = station_i->demand;

        // First verify if the change is feasible.
        if (qty_supplied_j + demand_i <= this->capacity) {

            this->InsertIntoRoute(
                    station_i,
                    station_j->route,
                    station_j);

            this->TwoOptExchange(station_i->route, CVRP::HARD_2_OPT);
            this->TwoOptExchange(station_j->route, CVRP::HARD_2_OPT);

            distance = this->total_distance;
        }
    // Operator (1,1), i.e. exchange between elements i and j.
    } else {

        if ((qty_supplied_i - demand_i) <= this->capacity
                && (qty_supplied_j - demand_j) <= this->capacity) {

            this->SwapStations(
                    station_i, 
                    route_i,
                    station_j,
                    route_j);

            this->TwoOptExchange(route_i, CVRP::HARD_2_OPT);
            this->TwoOptExchange(route_j, CVRP::HARD_2_OPT);

            distance = this->total_distance;
        }
    }

    //this->PrintRoutes();
    //this->Validate();

    return distance;
}

long CVRP::TwoOptExchange(CVRP::Route * route, CVRP::TwoOptMode mode) {

    long best_distance = route->distance;
    long new_distance = INT_MAX;

    /*printf("2-OPT: Starting local search on route: "\
            "\n\t[NUMBER] : %2d"\
            "\n\t[NS] : %3d"\
            "\n\t[D] : %3ld\n", 
            route->index,
            route->num_stations,
            route->distance);*/

    // For clearer code and readability, let's operate over an array, in 
    // which we could work with indexes.
    CVRP::Station** route_stations = (CVRP::Station**) calloc(route->num_stations + 2, sizeof(CVRP::Station*));
    CVRP::Station * curr_route_station = route->start_station;

    // Just a reference to the graph's depot.
    CVRP::Station * depot = this->GetDepot();

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

            /*printf("2-OPT: Evaluating potential of exchanging nodes "\
                    "(%2d [i = %d][%2d][%2d], %2d [j = %d][%2d][%2d]).\n", 
                    route_stations[i]->index, i, route_stations[i]->coords.x, route_stations[i]->coords.y,
                    route_stations[j]->index, j, route_stations[j]->coords.x, route_stations[j]->coords.y);*/

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

                /*printf("2-OPT: Produced new best distance by exchanging nodes "\
                        "(%2d [%d], %2d [%d]). \n\t[D] : %3ld\n", 
                        route_stations[i]->index, i, route_stations[j]->index, j,
                        new_distance);*/

                // Effectively swap the stations in the graph if the calling 
                // mode is different than CVRP::DRY_2_OPT.
                if (mode != CVRP::DRY_2_OPT) {
                    this->SwapStations(route_stations[i], route_stations[j], j - i);
                }

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
}

CVRP::IChange CVRP::RandomOneInterchange() {

    // Get pointers to 2 random stations (indexes between 1 and 
    // stations_size - 1), located in 2 different routes.
    int stations_size = this->stations.size();
    int i = 0;
    int j = 0;

    //long distance = 0;

    IChange ichange;

    // Keep trying until a change is accepted.
    bool ichange_accepted = false;

    do {

        do {

            // Never choose the depot, at stations[0].
            i = rand() % (stations_size - 1) + 1;
            j = rand() % (stations_size - 1) + 1;

        } while (stations[i]->route->index == stations[j]->route->index);

        // Choose between i-change operators (0) (0,1); (1) (1,0) or (2) 
        // (1,1).
        int opt = rand() % 3;
        //int opt = ICHANGE_OP11;

        int qty_supplied_i = stations[i]->route->qty_supplied;
        int qty_supplied_j = stations[j]->route->qty_supplied;

        int demand_i = stations[i]->demand;
        int demand_j = stations[j]->demand;

        // Pre-calculate the euclidean distances used for the 
        // route cost (distance) calculations of new solutions.
        int dist_i_j   = 
            CVRP::Euc2DistCalc(stations[i], stations[j]);
        int dist_pj_j  = 
            CVRP::Euc2DistCalc(stations[j]->prev, stations[j]);
        int dist_pi_i  = 
            CVRP::Euc2DistCalc(stations[i]->prev, stations[i]);
        int dist_pi_j  = 
            CVRP::Euc2DistCalc(stations[i]->prev, stations[j]);
        int dist_pj_nj = 
            CVRP::Euc2DistCalc(stations[j]->prev, stations[j]->next);
        int dist_j_ni  = 
            CVRP::Euc2DistCalc(stations[j], stations[i]->next);
        int dist_j_nj  = 
            CVRP::Euc2DistCalc(stations[j], stations[j]->next);
        int dist_pj_i  = 
            CVRP::Euc2DistCalc(stations[j]->prev, stations[i]);
        int dist_pi_ni = 
            CVRP::Euc2DistCalc(stations[i]->prev, stations[i]->next);
        int dist_i_ni  = 
            CVRP::Euc2DistCalc(stations[i], stations[i]->next);
        int dist_i_nj  = 
            CVRP::Euc2DistCalc(stations[i], stations[j]->next);

        // Operator (0,1), i.e. element j goes into route of element i.
        if (opt == ICHANGE_OP01) {

            // First verify if the change is feasible.
            if ((qty_supplied_i + demand_j) <= this->capacity) {

                ichange_accepted = true;

                ichange.distance = 
                    // Insert stations[j] after [i].
                    + dist_i_j
                    + dist_j_ni
                    - dist_i_ni

                    // New distance between stations[j] and [j + 1], 
                    // in route j.
                    + dist_pj_nj
                    - dist_pj_j
                    - dist_j_nj;

                ichange.operation = ICHANGE_OP01;
            }

        } else if (opt == ICHANGE_OP10) {

            // First verify if the change is feasible.
            if (qty_supplied_j + demand_i <= this->capacity) {

                ichange_accepted = true;

                ichange.distance = 
                    // Insert stations[i] after [j].
                    + dist_i_j
                    + dist_i_nj
                    - dist_j_nj

                    // New distance between stations[i] and [i + 1], 
                    // in route i.
                    + dist_pi_ni
                    - dist_pi_i
                    - dist_i_ni;

                ichange.operation = ICHANGE_OP10;
            }

        } else {

            if ((qty_supplied_i - demand_i + demand_j) <= this->capacity
                    && (qty_supplied_j - demand_j + demand_i) <= this->capacity) {

                ichange_accepted = true;

                ichange.distance = 
                    // Change in distance caused by including stations[i] in 
                    // route j, between stations [j - 1] and [j + 1].
                    dist_pj_i
                    + dist_i_nj
                    - dist_pi_i
                    - dist_i_ni

                    // Change in distance caused by including stations[j] in 
                    // route i, between stations [i - 1] and [i + 1].
                    + dist_pi_j
                    + dist_j_ni
                    - dist_pj_j
                    - dist_j_nj;

                ichange.operation = ICHANGE_OP11;
            }
        }

    } while (!ichange_accepted);

    ichange.i = i;
    ichange.j = j;
    ichange.distance += this->total_distance;

    // Perform the 1-interchange procedure.
    this->OneInterchange(ichange);

    return ichange;
}

int CVRP::GenerateRandomIndividual(
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

