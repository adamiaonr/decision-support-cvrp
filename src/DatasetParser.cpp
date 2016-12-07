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

/*
 * \brief   takes a line and breaks it into tokens
 * 
 * \param   line_to_tokenize        line to separate into tokens.
 * \param   delimiter               the delimiter between tokens (default is space ' ')
 * \param   tokens                  an array of strings, which will save the tokens.
 *
 * \return  the number of tokens resulting for the tokenization
 */
int tokenize_line(
    char * line_to_tokenize,
    const char * delimiter, 
    char ** tokens) {

    int n = 0;
    // to break the line into tokens, we use the strtok() function (see 
    // http://www.cplusplus.com/reference/cstring/strtok/). strtok() returns 
    // the first token in the string, and takes 2 arguments:
    //  -# the string to be broken into tokens (in our case, 'line')
    //  -# the character that separates each token (in our case, the space 
    //     character, ' ' or DELIMITER)
    //

    // we extract the tokens in 2 steps. let's use the line 
    // 'CAPACITY : 100' as an example.
    // step 1: we extract the 1st token of the line, i.e. 'CAPACITY'. this 
    // is done by calling strtok(line, DELIMITER)
    tokens[0] = strtok(line_to_tokenize, delimiter);

    // step 2: we extract the next tokens, i.e. ':' and '100'. to do 
    // so, we keep calling strtok(0, DELIMITER) until there are no more 
    // tokens, at which point strtok() returns NULL. 
    // if the line was empty (tokens[0] - the 1st token - is NULL), we don't 
    // call strtok() again
    if (tokens[0]) {

        // this is only here to avoid taking too many tokens per line
        for (n = 1; n < MAX_TOKENS_PER_LINE; n++) {

            // NOTE THAT AT STEP 2 WE DON'T CALL strtok() WITH 'line', WE USE 0 
            // INSTEAD. THIS IS BECAUSE strtok() REMEMBERS WE ARE BREAKING 'line'.
            // IF YOU WANTED strtok() TO START BREAKING A DIFFERENT STRING - SAY 'new_line'
            // YOU WOULD CALL strtok(new_line, DELIMITER).
            tokens[n] = strtok(0, delimiter);

            // if there are no more tokens, jump off the cycle
            if (!tokens[n]) 
                break;
        }
    }

    return n;
}

/*
 * \brief   opens a dataset file (a .vrp file) and extracts the CVRP information 
 *          into a CVRP object.
 *
 * we're interested in extracting the following information from the .vrp file:
 *  -# DIMENSION : number of nodes (stations) of the CVRP problem
 *  -# CAPACITY  : the capacity of the vehicle
 *  -# NODE_COORD_SECTION list : the <x,y> coordinates of the nodes
 *  -# DEMAND_SECTION list : the demands at each node
 * 
 * \param   dataset_filename        the name of the .vrp file to open and 
 *                                  extract the information from.
 * \param   graph                   the CVRP object (i.e. the problem graph) on 
 *                                  which the information will be saved. this 
 *                                  graph can then be used by an heuristic 
 *                                  to solve the problem (e.g. CWS).
 */
int DatasetParser::ParseDataset(
        const char * dataset_filename,
        CVRP * graph) {

    // carrying capacity for each vehicle
    int capacity = 0;
    // number of nodes (i.e. demand points or stations)
    int number_of_nodes = 0;
    // array of node descriptors, which will be used 
    DatasetParser::Node * nodes = NULL;

    // we use a std::ifstream object to read the dataset (.vrp) file
    ifstream fin;
    // open the file
    fin.open(dataset_filename);

    // if open() fails for some reason, abort
    if (!fin.good()) {

        fprintf(
                stderr, 
                "PARSE DATASET: Dataset file (%s) not found.\n",
                dataset_filename);

        return -1;
    }

    // read the file, line-by-line and repeat till the last line
    while (!fin.eof()) {

        // we read each line of the file into line
        char line[MAX_CHARS_PER_LINE];
        fin.getline(line, MAX_CHARS_PER_LINE);

        // the lines in the .vrp file are formatted in different ways.
        // 
        // e.g. to represent a simple parameter like capacity, the file uses 
        // a single like like:
        // "CAPACITY : 100"
        //
        // e.g. a list of arguments (the id and <x,y> coordinates of the 
        // stations) the file uses multiple lines, e.g.:
        // "NODE_COORD_SECTION 
        //  1 82 76 <- (i.e. id of station, x coordinate, y coordinate)
        //  2 96 44
        //  (...)""
        // 
        // note that for both cases, the different parts of each line are separated by spaces, i.e. 
        // the ' ' character (represented by DELIMITER, see line 43 in DataParser.h). 

        // we first break the line into smaller strings, i.e. 'tokens'. in the 
        // end, each 'token' will be an element in an array called 'tokens'. 
        // we define the function tokenize_line() to do so.
        char * tokens[MAX_TOKENS_PER_LINE] = {};
        int number_of_tokens = tokenize_line(line, DELIMITER, tokens);

        // now that we have the line divided in tokens, we look for the 
        // parameters that interest us and save them in variables to 
        // be used by our program. the variables of the .vrp file that 
        // are interesting to us are:
        //  -# DIMENSION
        //  -# CAPACITY
        //  -# NODE_COORD_SECTION list
        //  -# DEMAND_SECTION list
        for (int i = 0; i < number_of_tokens; i++) {

            // to determine the parameter referenced by the line, we look into 
            // the 1st token (tokens[0]), which contains the name of the parameter
            if (strcmp(tokens[i], "DIMENSION") == 0) {

                // we convert the 3rd token from a string to an integer (we 
                // skip the 2nd token, ':', which doesn't matter to us) 
                number_of_nodes = atoi(tokens[i + 2]);
                printf("DataParser::ParseDataset() : [INFO] number_of_nodes = %d\n", number_of_nodes);

            } else if (strcmp(tokens[i], "CAPACITY") == 0) {

                // same as above, now for CAPACITY
                capacity = atoi(tokens[i + 2]);
                printf("DataParser::ParseDataset() : [INFO] capacity = %d\n", capacity);

                // we update the 'capacity' attribute of the CVRP object
                graph->capacity = capacity;

            } else if (strcmp(tokens[i], "NODE_COORD_SECTION") == 0) {

                nodes = (DatasetParser::Node *) calloc (number_of_nodes, sizeof(DatasetParser::Node));

                // we've found the 'NODE_COORD_SECTION', below it we have 
                // n lines which have the <x,y> coordinates for each 
                // station. 

                // to retrieve the coordinates of station k (with k = {0, 1, 2, ...}), 
                // we read the lines which are below it, and tokenize each one 
                // the lines into 3 tokens: 'station id', 'x', 'y'
                for (int k = 0; k < number_of_nodes; k++) {

                    // read the line
                    fin.getline(line, MAX_CHARS_PER_LINE);
                    // tokenize the line into the coordinate_tokens
                    char * coordinate_tokens[MAX_TOKENS_PER_LINE] = {};
                    tokenize_line(line, DELIMITER, coordinate_tokens);
                    // update the coordinates for node k with the 2nd and 3rd 
                    // tokens (i.e. tokens[1] and tokens[2]: tokens[0] doesn't 
                    // matter here)
                    nodes[k].coords = { atoi(coordinate_tokens[1]), atoi(coordinate_tokens[2]) };
                    printf("DataParser::ParseDataset() : [INFO] nodes[%d].coords = { %d, %d }\n", 
                        k, nodes[k].coords.x, nodes[k].coords.y);
                }

            } else if (strcmp(tokens[i], "DEMAND_SECTION") == 0) {

                // same as above, but for the demands of each node
                for (int k = 0; k < number_of_nodes; k++) {

                    fin.getline(line, MAX_CHARS_PER_LINE);
                    // tokenize the line into the coordinate_tokens
                    char * demand_tokens[MAX_TOKENS_PER_LINE] = {};
                    tokenize_line(line, DELIMITER, demand_tokens);
                    // update the demand of node k with the 2nd token 
                    // (i.e. tokens[1]: tokens[0] doesn't matter here)
                    nodes[k].demand = atoi(demand_tokens[1]);
                    printf("DataParser::ParseDataset() : [INFO] nodes[%d].demand = %d\n", 
                        k, nodes[k].demand);
                }
            }
        }
    }

    // finally, initialize the CVRP graph, given the information retrieved from 
    // parsing the .vrp file, in the nodes array
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

    for (int k = 1; k < number_of_nodes; k++) {

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

