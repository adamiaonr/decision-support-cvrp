#include "DatasetParser.h"

//! DatasetParser constructor.
/*!
*/
DatasetParser::DatasetParser()
{
}

//! DatasetParser destructor.
/*!
*/
DatasetParser::~DatasetParser()
{
}

int DatasetParser::ParseDataset(
        const char * dataset,
        CVRP * graph) {

    int capacity = 0;
    int number_nodes = 0;

    DatasetParser::Node * nodes = NULL;

    ifstream fin;
    fin.open(dataset);

    if (!fin.good()) {

        fprintf(
                stderr, 
                "PARSE DATASET: Dataset file (%s) not found.\n",
                dataset);

        return -1;
    }

    while (!fin.eof()) {

        char buf[MAX_CHARS_PER_LINE];
        fin.getline(buf, MAX_CHARS_PER_LINE);

        // parse the line into blank-delimited tokens
        int n = 0; // a for-loop index

        // array to store memory addresses of the tokens in buf
        const char *token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0

        // parse the line
        token[0] = strtok(buf, DELIMITER); // first token

        if (token[0]) {

            for (n = 1; n < MAX_TOKENS_PER_LINE; n++) {

                token[n] = strtok(0, DELIMITER); // subsequent tokens

                if (!token[n]) break; // no more tokens
            }
        }

        // process (print) the tokens
        for (int i = 0; i < n; i++) {

            if (strcmp(token[i], "DIMENSION") == 0) {

                number_nodes = atoi(token[i + 2]);

                //printf("number_nodes = %d\n", number_nodes);

            } else if (strcmp(token[i], "CAPACITY") == 0) {

                capacity = atoi(token[i + 2]);

                //printf("capacity = %d\n", capacity);
                graph->capacity = capacity;

            } else if (strcmp(token[i], "NODE_COORD_SECTION") == 0) {

                nodes = (DatasetParser::Node *) calloc (number_nodes, sizeof(DatasetParser::Node));

                for (int k = 0; k < number_nodes; k++) {

                    fin.getline(buf, MAX_CHARS_PER_LINE);

                    int t = 0;

                    const char *n_token[MAX_TOKENS_PER_LINE] = {};
                    n_token[0] = strtok(buf, DELIMITER);

                    if (n_token[0]) {

                        for (t = 1; t < MAX_TOKENS_PER_LINE; t++) {

                            n_token[t] = strtok(0, DELIMITER);

                            if (!n_token[t]) break;
                        }
                    }

                    nodes[k].coords = { atoi(n_token[1]), atoi(n_token[2]) };

                    //printf("nodes[k].coords = { %d, %d }\n", nodes[k].coords.x, nodes[k].coords.y);
                }

            } else if (strcmp(token[i], "DEMAND_SECTION") == 0) {

                for (int k = 0; k < number_nodes; k++) {

                    fin.getline(buf, MAX_CHARS_PER_LINE);

                    int t = 0;

                    const char *n_token[MAX_TOKENS_PER_LINE] = {};
                    n_token[0] = strtok(buf, DELIMITER);

                    if (n_token[0]) {

                        for (t = 1; t < MAX_TOKENS_PER_LINE; t++) {

                            n_token[t] = strtok(0, DELIMITER);

                            if (!n_token[t]) break;
                        }
                    }

                    nodes[k].demand = atoi(n_token[1]);

                    //printf("nodes[k].demand = %d\n", nodes[k].demand);
                }
            }
        }
    }

    // Initialize the CVRP graph, given the information retrieved from 
    // parsing the .vrp file.

    if (graph->AddStation(
                1,
                nodes[0].demand,
                { nodes[0].coords.x, nodes[0].coords.y },
                true) < 0) {

        fprintf(
                stderr, 
                "ERROR: Depot addition to CVRP failed { %d, %d }.\n", 
                nodes[0].coords.x, 
                nodes[0].coords.y);
    }

    for (int k = 1; k < number_nodes; k++) {

        if (graph->AddStation(
                k + 1,
                nodes[k].demand,
                { nodes[k].coords.x, nodes[k].coords.y },
                false) < 0) {

            fprintf(
                stderr, 
                "ERROR: Station addition to CVRP failed { %d, %d }.\n", 
                nodes[k].coords.x, 
                nodes[k].coords.y);
        }
    }

    free(nodes);

    return 0;
}

