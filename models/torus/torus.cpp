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

    // for uniform random routing

//    this->registerRNG(this->rng_);
//    std::exponential_distribution<double> interval_expo(1.0 / this->mean_interval_);
//    std::uniform_int_distribution<unsigned int> random_direction(*this->rng); // random direction
//    auto packet_arrival = (direction_t) random_direction(*this->rng);
//
//    // this doesn't work, needs to be configurable for dimensionality
//    std::uniform_int_distribution<int> random_x(-99,100);
//    std::uniform_int_distribution<int> random_y(-99,100);
//    std::uniform_int_distribution<int> random_z(-99,100);
//
//    for (unsigned int i = 0; i < this->total_nodes; i++) {
//        // new MessageEvent(name, arrival, random_x(rng_), random_y(rng_), packet_arrival, (unsigned int) std::ceil(interval_expo(*this->rng_)));
//        events.emplace_back(new MessageEvent {});
//    }

    return events;
}

std::vector<std::shared_ptr<warped::Event> > Node::receiveEvent(const warped::Event& event) {
    std::vector<std::shared_ptr<warped::Event>> events;
    auto message_event = static_cast<const MessageEvent&>(event);
    unsigned int event_ts = message_event.event_ts_ + TS_INTERVAL;
    message_event.hop_count_++;

    // check that a message has reached it's destination
    if (message_event.receiver_name_ == node_name(index_)) {
//        events.emplace_back(new MessageEvent { node_name(index_), event_ts });
    } else {
        // else, forward on
//        events.emplace_back(new MessageEvent { node_name(neighbor(index_, message_event.destination)), event_ts });
    }

    // wait queue?

    switch (message_event.type_) {
      case GENERATE: {
          // state, bf - buffer?, message, lp
          int destination = message_event.destination_index;
          int grid_dimension = grid_dimension_;

          if(destination < grid_dimension) {
            message_event.destination_index = state_.neighbor_minus_index[destination];
          } else if(destination >= grid_dimension && destination < 2 * grid_dimension) {
            message_event.destination_index = state_.neighbor_minus_index[destination - grid_dimension];
          }
        } break;
      case ARRIVAL: {

          // traffic
          // switch on current lane
              // increment counter in state
              // set [event type], emplace back with new variables, [event type], new x and y, timestamp, arrival from

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

// unsigned int Node::route(MessageEvent *message) {
unsigned int Node::route(unsigned int destination_index, unsigned int source_index, unsigned int direction) {
    // direction depends on the order of the torus
    // check the shortest path based on destination
    // dimension order routing
    // return a new index
    // get node from indices
    Node destination;  // NPE
    Node source;  // NPE
    int dimensions[grid_dimension];
    int dest[grid_dimension];
    dimensions[0] = *destination;

    // find destination dimensions using destination LP ID
    for (int i = 0; i < grid_dimension; i++ ) {
        dest[i] = dimensions[i] % dim_length[i];
        dimensions[i + 1] = (dimensions[i] - dest[i]) / dim_length[i];
    }

    for(int i = 0; i < grid_dimension; i++ ) {
        if (state_->dim_position[i] - dest[i] > half_length[i])  {
            message->destination_index = state_->neighbor_minus_index[i];
            message->destination_dimension = i;
            message->destination_direction = 1;
        } else if (state_->dim_position[i] - dest[i] < -half_length[i]) {
            message->destination_index = state_->neighbor_minus_index[i];
            message->destination_dimension = i;
            message->destination_direction = 0;
        } else if ((state_->dim_position[i] - dest[i] <= half_length[i]) && (state_->dim_position[i] - dest[i] > 0 )) {
            message->destination_index = state_->neighbor_minus_index[i];
            message->destination_dimension = i;
            message->destination_direction = 0;
        } else if ((state_->dim_position[i] - dest[i] >= -half_length[i]) && (state_->dim_position[i] - dest[i] < 0)) {
            message->destination_index = state_->neighbor_minus_index[i];
            message->destination_dimension = i;
            message->destination_direction = 1;
        }
    }

    return 1;
}

int main(int argc, const char **argv) {
    unsigned int buffer_slots = 1000;

    /* Set the default values for arguments */
    unsigned int grid_dimension = 5;
    unsigned int grid_size = 1000;  // configurable per dimension in ROSS
    unsigned int grid_order = 4;
    unsigned int mean_interval = 200;
    // ROSS sets opt_memory, mpi_message_size, mem_factor, num_mpi_msgs
    // add parameter for different routing algorithms
    /* Read arguments */
    TCLAP::ValueArg<unsigned int> grid_dimension_arg("d", "dimension",
                    "Dimensionality of the torus", false, grid_dimension, "unsigned int");
    TCLAP::ValueArg<unsigned int> grid_size_arg("s", "size",
                    "Size of the torus grid", false, grid_size, "unsigned int");
    TCLAP::ValueArg<unsigned int> grid_order_arg("o", "order",
                    "Order of nodes in torus", false, grid_order, "unsigned int");
    TCLAP::ValueArg<unsigned int> mean_interval_arg("i", "mean-interval",
                    "Mean interval", false, mean_interval, "unsigned int");

    std::vector<TCLAP::Arg*> cmd_line_args = { &grid_dimension_arg,
                                               &grid_size_arg,
                                               &grid_order_arg,
                                               &mean_interval_arg
                                             };

    warped::Simulation simulation {"Torus Network Simulation", argc, argv, cmd_line_args};

    grid_dimension = grid_dimension_arg.getValue();
    grid_size      = grid_size_arg.getValue();
    grid_order     = grid_order_arg.getValue();
    mean_interval  = mean_interval_arg.getValue();

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
    // nlp_nodes_per_pe = total_nodes/tw_nnodes()/g_tw_npe; ??
    unsigned int rows = sqrt(total_nodes);  // need to track total nodes
    unsigned int cols = rows;
    // total_lps
    // node_rem - ROSS mapping
    unsigned int num_packets = 1;
    unsigned int num_chunks = PACKET_SIZE / chunk_size;
//    g_tw_mapping=CUSTOM;
//    g_tw_custom_initial_mapping=&torus_mapping;
//    g_tw_custom_lp_global_to_local_map=&torus_mapping_to_lp;

    // create dimensions array, and initialize to each LP
    unsigned int dimension_coordinates[grid_dimension];  // +1

    /* Torus is a grid of size n with k dimensions */
    for (unsigned int k = 0; k < grid_dimension; k++) {
        for (unsigned int n = 0; n < grid_size; n++) {
            lps.emplace_back(
                    node_name(n),
                    grid_dimension,
                    grid_size,
                    grid_order,
                    mean_interval,
                    k,
                    n
            );
        }
    }

    for (auto& lp : lps) {
        lp_pointers.push_back(&lp);
    }

    // have to separate initialization steps
    for (unsigned int i = 0; i < grid_dimension; i++) {
        // initialize each to lp->gid ??

        // find each LP's coordinates

        // real coordinates (ignoring sim array)
        // node_state.position[i] = dim[i] % length;
        // dim[i+1] = i - node_state.position / length;
        // half_length = length / 2;

        // initialize "factor" to calculate neighbors
        unsigned int factor[];
        // factor [dimension[ = 1, for j<dimension, factor *= size;

        // calculate dimension neighbors
        unsigned int dimension_plus_coordinates[];
        unsigned int dimension_minus_coordinates[];

        // set coordinates[i] = node_state.position[i]

        // calculate +/- 1 neighbor's LP index

        //        temp_dim_plus_pos[j] = (state_->dim_position[j] + 1 + dim_length[j]) % dim_length[j];
//        temp_dim_minus_pos[j] =  (state_->dim_position[j] - 1 + dim_length[j]) % dim_length[j];
//
//        state_->neighbour_minus_lpID[j] = 0;
//        state_->neighbour_plus_lpID[j] = 0;
//
//        for ( i = 0; i < grid_dimension; i++ )
//        {
//            state_->neighbour_minus_lpID[j] += factor[i] * temp_dim_minus_pos[i];
//            state_->neighbour_plus_lpID[j] += factor[i] * temp_dim_plus_pos[i];
//        }
//
//        temp_dim_plus_pos[j] = state_->dim_position[j];
//        temp_dim_minus_pos[j] = state_->dim_position[j];
    }

// initialize lps
//    for( j=0; j < 2 * grid_dimension; j++ )
//    {
//        for( i = 0; i < NUM_VC; i++ )
//        {
//            state_->buffer[j][i] = 0;
//            state_->next_link_available_time[j][i] = 0.0;
//        }
//    }
//    // record LP time
//    state_->packet_counter = 0;
//    state_->waiting_list = tw_calloc(TW_LOC, "waiting list", sizeof(struct waiting_packet), WAITING_PACK_COUNT);
//
//    for (j = 0; j < WAITING_PACK_COUNT - 1; j++) {
//        state_->waiting_list[j].next = &state_->waiting_list[j + 1];
//        state_->waiting_list[j].dim = -1;
//        state_->waiting_list[j].dir = -1;
//        state_->waiting_list[j].packet = NULL;
//    }
//
//    state_->waiting_list[j].next = NULL;
//    state_->waiting_list[j].dim = -1;
//    state_->waiting_list[j].dir = -1;
//    state_->waiting_list[j].packet = NULL;
//
//    state_->head = &state_->waiting_list[0];
//    state_->wait_count = 0;

    simulation.simulate(lp_pointers);

    return 0;
}
