#include <vector>
#include <memory>
#include <cassert>
#include <algorithm>
#include <random>
#include <cstdlib>

#include "warped.hpp"
#include "torus.hpp"
#include "ppm/ppm.hpp"
#include "tclap/ValueArg.h"

#define TS_INTERVAL 1

WARPED_REGISTER_POLYMORPHIC_SERIALIZABLE_CLASS(MessageEvent)

std::string node_name (unsigned int index) {
    return std::string("Node_") + std::to_string(index);
}

std::vector<std::shared_ptr<warped::Event> > Node::initializeLP() {
    std::vector<std::shared_ptr<warped::Event>> events;

    events.emplace_back( new MessageEvent { node_name(thistate->index_), TS_INTERVAL } );

    return events;
}

std::vector<std::shared_ptr<warped::Event> > Node::receiveEvent(const warped::Event& event) {
    std::vector<std::shared_ptr<warped::Event>> events;
    auto message_event = static_cast<const MessageEvent&>(event);
    unsigned int event_ts = message_event.event_ts_ + TS_INTERVAL;

    message_event.hop_count_++;

    // check that a message has reached it's destination
    if (message_event.receiver_name_ == node_name(thistate->index_)) {
        events.emplace_back(new MessageEvent { node_name(thistate->index_), event_ts });
    } else {
        // else, forward on
        events.emplace_back(new MessageEvent { node_name(neighbor(thistate->index_, message_event.destination)), event_ts });
    }

    // while there is an event available?
    switch (thistate->type_) {
      case GENERATE: {

        } break;
      case ARRIVAL: {

        } break;
      case SEND: {

        } break;
      case PROCESS: {

        } break;
      case CREDIT: {

        } break;
      case WAIT: {

        } break;
      case MPI_SEND: {

        } break;
      case MPI_RECV: {

        } break;
    }


    return events;
}

unsigned int Node::neighbor(unsigned int destination_index,
                            unsigned int source_index,
                            unsigned int direction) {
    // direction depends on the order of the torus

    // check the shortest path based on destination

    // dimension order routing

    // return a new index
    // return 1;

    // get node from indices
    Node destiation;
    Node source;

    int dim_N[ grid_dimension ],
    dest[ grid_dimension ],
    i;

    dim_N[0] = *destination;

    // find destination dimensions using destination LP ID
    for ( i = 0; i < grid_dimension; i++ ) {
        dest[ i ] = dim_N[ i ] % dim_length[ i ];
        dim_N[ i + 1 ] = ( dim_N[ i ] - dest[ i ] ) / dim_length[ i ];
    }

    for( i = 0; i < grid_dimension; i++ ) {
        if ( state->dim_position[ i ] - dest[ i ] > half_length[ i ] )  {
            *destination = state->neighbour_plus_lpID[ i ];
            *dim = i;
            *direction = 1;
            break;
        }

        if ( state->dim_position[ i ] - dest[ i ] < -half_length[ i ] ) {
            *destination = state->neighbour_minus_lpID[ i ];
            *dim = i;
            *direction = 0;
            break;
        }

        if ( ( state->dim_position[ i ] - dest[ i ] <= half_length[ i ] ) && ( state->dim_position[ i ] - dest[ i ] > 0 ) ) {
            *destination = state->neighbour_minus_lpID[ i ];
            *dim = i;
            *direction = 0;
            break;
        }

        if (( state->dim_position[ i ] - dest[ i ] >= -half_length[ i ] ) && ( state->dim_position[ i ] - dest[ i ] < 0) ) {
            *destination = state->neighbour_plus_lpID[ i ];
            *dim = i;
            *direction = 1;
            break;
        }
    }
}

int main(int argc, const char **argv) {
    unsigned int buffer_slots = 1000;

    /* Set the default values for arguments */
    unsigned int grid_dimension = 5;
    unsigned int grid_size = 1000;  // configurable per dimension in ROSS
    unsigned int grid_order = 4;

    // ROSS sets opt_memory, mpi_message_size, mem_factor, num_mpi_msgs

    // add parameter for different routing algorithms

    /* Read arguments */
    TCLAP::ValueArg<unsigned int> grid_dimension_arg("d", "dimension",
                    "Dimensionality of the torus", false, grid_dimension, "unsigned int");
    TCLAP::ValueArg<unsigned int> grid_size_arg("s", "size",
                    "Size of the torus grid", false, grid_size, "unsigned int");
    TCLAP::ValueArg<unsigned int> grid_order_arg("o", "order",
                    "Order of nodes in torus", false, grid_order, "unsigned int");

    std::vector<TCLAP::Arg*> cmd_line_args = {   &grid_dimension_arg,
                                                 &grid_size_arg,
                                                 &grid_order_arg
                                             };

    warped::Simulation simulation {"Torus Network Simulation", argc, argv, cmd_line_args};

    grid_dimension = grid_dimension_arg.getValue();
    grid_size = grid_size_arg.getValue();
    grid_order = grid_order_arg.getValue();

    unsigned int MEAN_INTERVAL = 200;
    unsigned int injection_limit = 10;
    unsigned int injection_interval = 20000;
    unsigned int vc_size = 16384;
    unsigned int link_bandwidth = 2;

    // ross global variables
    unsigned int PACKET_SIZE = 512; // why packet size?
    unsigned int chunk_size = 32;
    unsigned int NUM_VC = 1; // constant for VC? virtual channel?
    unsigned int COLLECT_POINTS = 100; // number of points to collect stats from
    unsigned int WAITING_PACK_COUNT = 1; // buffers can be full

    // tw_stime tracking simulation time statistics
    unsigned int N_finished_packets = 0;
    unsigned int N_finished_msgs = 0;
    /* number of finished packets, generated packets, size of buffer queues
     * and average number of hops travelled by the packets at different
     * points during the simulation per PE */
    unsigned int N_finished_storage[COLLECT_POINTS];
    unsigned int N_generated_storage[COLLECT_POINTS];
    unsigned int N_queue_depth[COLLECT_POINTS];
    unsigned int N_num_hops[COLLECT_POINTS];
    unsigned int total_hops = 0;
    unsigned int half_length[grid_dimension_];  // ROSS - for calculating dimensions of coordinates for a node
    unsigned int half_length_sim[grid_dimension_sim];
    // ROSS simulation statistics removed

    /* Create the LPs */
    std::vector<Node> lps;
    std::vector<warped::LogicalProcess*> lp_pointers;

    unsigned int nodes = grid_dimension * grid_size;
    // nlp_nodes_per_pe = N_nodes/tw_nnodes()/g_tw_npe; ??
    unsigned int rows = sqrt(N_nodes);
    unsigned int cols = rows;
    // total lps
    // node_rem
    num_packets = 1;
    num_chunks = PACKET_SIZE / chunk_size;
//    g_tw_mapping=CUSTOM;
//    g_tw_custom_initial_mapping=&torus_mapping;
//    g_tw_custom_lp_global_to_local_map=&torus_mapping_to_lp;

    // create dimensions array, and initialize to each LP
    unsigned int dimension_coordinates[grid_dimension];  // +1

    // have to separate initialization steps
    for (int i = 0; i < grid_dimension; i++) {
        // initialize each to lp->gid ??

        // find each LP's coordinates

        // real coordinates (ignoring sim array)
        // node_state.position[i] = dim[i] % length;
        // dim[i+1] = i - node_state.position / length;
        // half_length = length / 2;

        // initialize "factor" to calculate neighbors
        unsigned int factor[] = 1;
        // factor [dimension[ = 1, for j<dimension, factor *= size;

        // calculate dimension neighbors
        unsigned int dimension_plus_coordinates[];
        unsigned int dimension_minus_coordinates[];

        // set coordinates[i] = node_state.position[i]

        // calculate +/- 1 neighbor's LP index

//        temp_dim_plus_pos[ j ] = (state->dim_position[ j ] + 1 + dim_length[ j ]) % dim_length[ j ];
//        temp_dim_minus_pos[ j ] =  (state->dim_position[ j ] - 1 + dim_length[ j ]) % dim_length[ j ];
//
//        state->neighbour_minus_lpID[ j ] = 0;
//        state->neighbour_plus_lpID[ j ] = 0;
//
//        for ( i = 0; i < N_dims; i++ )
//        {
//            state->neighbour_minus_lpID[ j ] += factor[ i ] * temp_dim_minus_pos[ i ];
//            state->neighbour_plus_lpID[ j ] += factor[ i ] * temp_dim_plus_pos[ i ];
//        }
//
//        temp_dim_plus_pos[ j ] = state->dim_position[ j ];
//        temp_dim_minus_pos[ j ] = state->dim_position[ j ];
    }

    // initialize lps
//    for( j=0; j < 2 * N_dims; j++ )
//    {
//        for( i = 0; i < NUM_VC; i++ )
//        {
//            state->buffer[ j ][ i ] = 0;
//            state->next_link_available_time[ j ][ i ] = 0.0;
//        }
//    }
//    // record LP time
//    state->packet_counter = 0;
//    state->waiting_list = tw_calloc(TW_LOC, "waiting list", sizeof(struct waiting_packet), WAITING_PACK_COUNT);
//
//    for (j = 0; j < WAITING_PACK_COUNT - 1; j++) {
//        state->waiting_list[j].next = &state->waiting_list[j + 1];
//        state->waiting_list[j].dim = -1;
//        state->waiting_list[j].dir = -1;
//        state->waiting_list[j].packet = NULL;
//    }
//
//    state->waiting_list[j].next = NULL;
//    state->waiting_list[j].dim = -1;
//    state->waiting_list[j].dir = -1;
//    state->waiting_list[j].packet = NULL;
//
//    state->head = &state->waiting_list[0];
//    state->wait_count = 0;

    /* Torus is a grid of size n with k dimensions */
//    for (unsigned int k = 0; k < grid_dimension; k++) {
//        for (unsigned int n = 0; n < grid_size; n++) {
//            lps.emplace_back(node_name(n), grid_dimension, grid_size, grid_order, k, n);
//        }
//    }
//
//    for (auto& lp : lps) {
//        lp_pointers.push_back(&lp);
//    }

    simulation.simulate(lp_pointers);

    return 0;
}
