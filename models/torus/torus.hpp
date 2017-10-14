#ifndef TORUS_HPP
#define TORUS_HPP

#include <string>
#include <memory>
#include <random>

#include "warped.hpp"

// traffic enum for multiple routing algorithms
enum traffic {
    //UNIFORM_RANDOM=1,
    NEAREST_NEIGHBOR,
    //DIAGONAL
};

// mpi process struct has row and column for each mpi process for matrix transpose traffic
// struct mpi_process
// {
//  unsigned long long message_counter;
//  tw_stime available_time;
//
//  /*For matrix transpose traffic, we have a row and col value for each MPI process so that the message can be sent to the corresponding transpose of
//  the MPI process */
//  int row, col;
// };

// event types
enum message_event_t {
    GENERATE = 1,
    ARRIVAL,
    SEND,
    PROCESS,
    CREDIT,
    WAIT,
    MPI_SEND,
    MPI_RECV
};

WARPED_DEFINE_LP_STATE_STRUCT(NodeState) {
    unsigned int packet_counter;
    unsigned int next_link_available_ts[][]; // Time when the VC will be available for sending packets on this node
    unsigned int next_credit_available_ts[][]; // Time when the next VC will be available for sending credits on this node
    unsigned int buffer[][]; // Buffer occupancy of the current VC
    unsigned int dim_position[]; // torus dimension coordinates of this node
    unsigned int neighbor_minus_index[]; // torus neighbor coordinates for this node
    unsigned int neighbor_plus_index[]; // torus plus neighbor coordinates for this node, **removed simulated neighbor coordinates
    unsigned int source_dim;
    unsigned int direction;
    struct waiting_packet *waiting_list; // first element of linked list
    struct waiting_packet *ead;
    unsigned int wait_count;
};

// Linked list for storing waiting packets in the queue
struct waiting_packet {
    unsigned int dim;
    unsigned int dir;
    MessageEvent *packet;
    struct waiting_packet *next;
};

class MessageEvent : public warped::Event {
public:
    MessageEvent() = default;
    MessageEvent( const std::string receiver_name,
                  const message_event_t type,
                  unsigned int event_ts  )
        :   receiver_name_(receiver_name),
            type_(type),
            event_ts_(event_ts) {}

    const std::string& receiverName() const { return receiver_name_; }
    unsigned int timestamp() const { return event_ts_; }
    unsigned int size() const {
        return receiver_name_.length() + sizeof(event_ts_);
    }

    std::string receiver_name_;
    message_event_t type_;
    unsigned int hop_count_; // my_N_hop
    unsigned int event_ts_;

    // unsigned int chunk_id; // chunk ID of the packet
    // event time: mpi_send, mpi_recv, packet_generate etc.
    unsigned int source_dimension; // originating torus node dimension and direction
    unsigned int source_direction;
    unsigned int source_index;
    unsigned int destination[grid_dimension_]; // destination torus coordinates for this packet/message
    unsigned int destination_dimension;
    unsigned int destination_direction;
    unsigned int destination_index;
    unsigned int next_stop; // tw_lpid next stop of this message/packet
    unsigned int packet_ID;
    unsigned int packet_size; // size in bytes
    unsigned int wait_dir; // For waiting messages/packets that don't get a slot in the buffer
    unsigned int wait_dim;
    unsigned int wait_loc; // for packets waiting to be injected into the network
    unsigned int wait_type;
    unsigned int travel_start_time; // time when the packet starts travelling

    WARPED_REGISTER_SERIALIZABLE_MEMBERS(cereal::base_class<warped::Event>(this),
                                            receiver_name_, event_ts_)
};


class Node : public warped::LogicalProcess {
public:
    Node( const std::string& name,
          const unsigned int mean_interval,
          const unsigned int grid_dimension,
          const unsigned int grid_size,
          const unsigned int grid_order,
          const unsigned int dimension_index,
          const unsigned int index )
      : LogicalProcess(name),
        state_(),
        mean_interval_(mean_interval),
        grid_dimension_(grid_dimension),
        grid_size_(grid_size),
        grid_order_(grid_order),
        dimension_index_(dimension_index),
        index_(index) {

        state_.packet_counter;
        state_.next_link_available_ts = new unsigned int[2 * grid_dimension_][NUM_VC];
        state_.next_credit_available_ts = new unsigned int[2 * grid_dimension_][NUM_VC];
        state_.buffer = new unsigned int[2 * grid_dimension_][NUM_VC];
        state_.dim_position = new unsigned int[grid_dimension_];
        state_.neighbor_minus_index = new unsigned int[grid_dimension_];
        state_.neighbor_plus_index = new unsigned int[grid_dimension_];
        state_.source_dim;
        state_.direction;
        state_.waiting_list; // first element of linked list
        state_.read;
        state_.wait_count;
    }

    virtual warped::LPState& getState() { return state_; }
    virtual std::vector<std::shared_ptr<warped::Event>> initializeLP() override;
    virtual std::vector<std::shared_ptr<warped::Event>> receiveEvent(const warped::Event&);

    NodeState state_;

    // N_dims, removed simulation dimensions
//    unsigned int num_packets; // ROSS globals
//    unsigned int num_chunks;
    unsigned int packet_offset = 0;
    unsigned int num_buf_slots;
    unsigned int factor[grid_dimension_]; // for calculating torus dimensions // unsigned int factor_sim[grid_dimension_sim];
    float head_delay = 0.0; // calculating delays using the link bandwidth
    float credit_delay = 0.0;

protected:
    const unsigned int mean_interval_;
    const unsigned int grid_dimension_;
    const unsigned int grid_size_;
    const unsigned int grid_order_;
    const unsigned int dimension_index_;
    const unsigned int index_;
    unsigned int route(unsigned int destination_index, unsigned int source_index, unsigned int direction);
};

#endif
