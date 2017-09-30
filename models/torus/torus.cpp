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

    events.emplace_back( new MessageEvent { node_name(this->index_), TS_INTERVAL } );

    return events;
}

std::vector<std::shared_ptr<warped::Event> > Node::receiveEvent(const warped::Event& event) {
    std::vector<std::shared_ptr<warped::Event>> events;
    auto message_event = static_cast<const MessageEvent&>(event);
    unsigned int event_ts = message_event.event_ts_ + TS_INTERVAL;

    // ROSS credits the sender

    message_event.hop_count_++;

    // check that a message has reached it's destination
    if (message_event.receiver_name_ == node_name(this->index_)) {
        events.emplace_back(new MessageEvent { node_name(this->index_), event_ts });
    } else {
        // else, forward on
        events.emplace_back(new MessageEvent { node_name(neighbor(this->index_, message_event.destination)), event_ts });
    }

    return events;
}

unsigned int Node::neighbor(unsigned int index, unsigned int destination) {
    // direction depends on the order of the torus

    // check the shortest path based on destination

    index++;
    destination++;

    // return a new index
    return 1;
}

int main(int argc, const char **argv) {
    /* Set the default values for arguments */
    unsigned int grid_dimension = 5;
    unsigned int grid_size = 1000;
    unsigned int grid_order = 4;

    // run time arguments
    static int opt_mem = 3000;
    static long mpi_message_size = 32;
    static int num_mpi_msgs = 50;
    static int mem_factor = 16;
    static double MEAN_INTERVAL=200.0;
    static int TRAFFIC = UNIFORM_RANDOM;
    static int injection_limit = 10;
    static double injection_interval = 20000;
    static int vc_size = 16384;
    static double link_bandwidth = 2.0;
    static char traffic_str[512] = "uniform";

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

    /* Create the LPs */
    std::vector<Node> lps;
    std::vector<warped::LogicalProcess*> lp_pointers;

    /* Torus is a grid of size n with k dimensions */
    for (unsigned int k = 0; k < grid_dimension; k++) {
        for (unsigned int n = 0; n < grid_size; n++) {
            lps.emplace_back(node_name(n), grid_dimension, grid_size, grid_order, k, n);
        }
    }

    for (auto& lp : lps) {
        lp_pointers.push_back(&lp);
    }

    simulation.simulate(lp_pointers);

    return 0;
}
