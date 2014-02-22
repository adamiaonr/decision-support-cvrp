/*! \file   CVRP.h
    \brief  CVRP header file.
*/

#ifndef CVRP_H_
#define CVRP_H_

#include <list>
#include <cmath>
#include <vector>
#include <string>
#include <bitset>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <algorithm>
#include <sys/types.h>
//#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace std;
//using namespace boost;

/*! \def    ICHANGE_OP01
    \brief  
*/
#define ICHANGE_OP01 0

/*! \def    ICHANGE_OP10
    \brief  
*/
#define ICHANGE_OP10 1

/*! \def    ICHANGE_OP11
    \brief  
*/
#define ICHANGE_OP11 2

class CVRP {

    public:

        CVRP();
        CVRP(std::string dataset);
        ~CVRP();

        struct Station;
        struct Route;

        std::string dataset;

        uint32_t id;

        int num_stations;
        int num_routes;
        int capacity;

        long total_distance;

        Station * depot;

        // Choosing vector is safe, since without station removals there's no 
        // danger of iterator invalidations.
        std::vector <Station *> stations;

        // Must be list, due to route removal and possible iterator 
        // invalidations.
        std::list <Route *> routes;

        enum StationType {NORMAL, DEPOT};

        typedef struct euc_2d_coords { 

            int x;
            int y;

        } euc_2d_coords;

        typedef struct Route { 

            int index;
            int qty_supplied;
            int capacity;
            int num_stations;

            long distance;

            bool active;

            Station * start_station;
            Station * end_station;

        } Route;

        typedef struct Station { 

            // Index of the station, as defined in the CVRP datasets (typically 
            // from 1 - depot - to N).
            int index;

            // Demand qty for this particular station.
            int demand;

            // Euclidean 2D coordinates for this station.
            euc_2d_coords coords;

            // Type can be NORMAL or DEPOT.
            StationType type;

            // Pointers to next and previous (adjacent) stations.
            Station * next;
            Station * prev;

            // Pointer to the route this station belongs to.
            Route * route;

        } Station;

        typedef struct IChange {

            int i;
            int j;

            long distance;

            int operation;

        } IChange;

        static long Euc2DistCalc(euc_2d_coords i, euc_2d_coords j);
        static long Euc2DistCalc(Station * i, Station * j);

        int FillFromDataset(const char * dataset);

        void CopyFrom(CVRP * graph);

        Station * AddStation(
                int index,
                int demand,
                euc_2d_coords coords,
                bool is_depot);

        Station * GetNextStation(Station * station);

        Station * GetPrevStation(Station * station);

        void LinkStations(Station * start, Station * end);

        long SwapStations(Station * i, Station * j, int num_stations);

        long SwapStations(
                Station * i,
                Route * route_i,
                Station * j,
                Route * route_j);

        //void DeleteStation(Station * station);

        void PrintStations();

        Station * GetDepot();

        Route * AddRoute(
                int qty_supplied,
                long distance,
                Station * start_station,
                Station * end_station);

        Station** CreateRouteRep(Route * route);

        void DeleteRoute(Route * route);

        void PrintRoutes();

        Route * Assign2Route(Station * station, Route * route);

        long InsertIntoRoute(
                CVRP::Station * i,
                CVRP::Route * route_j,
                CVRP::Station * before_j);

        enum TwoOptMode {DRY_2_OPT, HARD_2_OPT};

        long TwoOptExchange(CVRP::Route * route, CVRP::TwoOptMode mode);

        long OneInterchange(CVRP::IChange & ichange);
        CVRP::IChange RandomOneInterchange();

        bool Validate();
        long CalcRouteDist(CVRP::Route * rt);

        static int GenerateRandomIndividual(
                CVRP * & individual,
                std::string dataset, 
                int num_routes);
};

#endif  // CVRP_H_
