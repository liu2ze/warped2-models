#ifndef TORUS_HPP
#define TORUS_HPP

#include <string>
#include <memory>
#include <random>

#include "warped.hpp"

using namespace std;

unsigned long N_dims = 5;
unsigned long N_dims_sim = 5;
unsigned long PACKET_SIZE = 512;
unsigned long NUM_VC = 1;
unsigned long N_COLLECT_POINTS = 100;;
static int dim_length_sim[] = {4,4,4,4,2};// 512 node case

// waiting packet count -- buffers can be full
// traffic enum for multiple routing algorithms
// enum traffic
// {
//   UNIFORM_RANDOM=1,
//   NEAREST_NEIGHBOR,
//   DIAGONAL
// };

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
    bool available_ = true;
    unsigned long long packet_counter;
    int buffer[2*N_dims][NUM_VC]; // Buffer occupancy of the current VC
    int dim_position[N_dims]; // torus dimension coordinates of this nod
    int dim_position_sim[N_dims_sim]; // torus dimension coordinates of the simulated torus dimension by this node (For TOPC paper)
    int neighbour_minus_lpID[N_dims]; // torus neighbor coordinates for this node
    int neighbour_plus_lpID[N_dims]; // torus plus neighbor coordinates for this node

    // For simulation purposes: Making the same nearest neighbor traffic
    // across all torus dimensions
    /* neighbors of the simulated torus dimension for this node */
    int neighbour_minus_lpID_sim[N_dims_sim];
    int neighbour_plus_lpID_sim[N_dims_sim];
    int source_dim;
    int direction;
    //first element of linked list
    struct waiting_packet *waiting_list;
    struct waiting_packet *ead;
    long wait_count;
};

// Linked list for storing waiting packets in the queue
struct waiting_packet {
    int dim;
    int dir;
    nodes_message *packet;
    struct waiting_packet *next;
};

// tw_stime tracking simulation time statistics

unsigned long N_finished_packets = 0;
unsigned long N_finished_msgs = 0;
/* number of finished packets, generated packets, size of buffer queues
 * and average number of hops travelled by the packets at different
 * points during the simulation per PE */
unsigned long N_finished_storage[N_COLLECT_POINTS];
unsigned long N_generated_storage[N_COLLECT_POINTS];
unsigned long N_queue_depth[N_COLLECT_POINTS];
unsigned long	N_num_hops[N_COLLECT_POINTS];
unsigned long total_hops = 0;
/* for calculating torus dimensions of real and simulated torus coordinates of
 * a node*/
int half_length[N_dims];
int half_length_sim[N_dims_sim];
// ROSS simulation statistics removed
int num_packets;/* number of packets in a message and number of chunks in a packet */
int num_chunks;
int packet_offset = 0;
const int chunk_size = 32;
int num_buf_slots;
int node_rem = 0;  // for ROSS mapping purposes
int num_rows, num_cols;
int factor[N_dims];/* for calculating torus dimensions */
int factor_sim[N_dims_sim];
float head_delay=0.0;/* calculating delays using the link bandwidth */
float credit_delay = 0.0;

class MessageEvent : public warped::Event {
public:
    MessageEvent() = default;

    MessageEvent( const string   receiver_name,
                  unsigned int        event_ts  )

        :   receiver_name_(receiver_name),
            event_ts_(event_ts) {}

    const string& receiverName() const { return receiver_name_; }

    unsigned int timestamp() const { return event_ts_; }

    unsigned int size() const {
        return receiver_name_.length() + sizeof(event_ts_);
    }

    string receiver_name_;
    string destination_name;
    unsigned int destination;
    unsigned int hop_count_;
    unsigned int event_ts_;

    unsigned long long packet_ID; // packet ID of the packet
    nodes_event_t type; // event time: mpi_send, mpi_recv, packet_generate etc.
    unsigned int source_dim; // originating torus node dimension and direction
    unsigned int source_direction;
    unsigned int saved_src_dim; // originating torus node dimension and direction (for reverse computation)
    unsigned int saved_src_dir;
    unsigned int wait_dir; // For waiting messages/packets that don't get a slot in the buffer
    unsigned int wait_dim;
    unsigned int dest[N_dims]; // destination torus coordinates for this packet/message
    unsigned int my_N_hop; // number of hops travelled by this message
    tw_lpid next_stop; // next stop of this message/packet
    unsigned int packet_size; // size of the message/packet in bytes
    short chunk_id; // chunk ID of the packet
    unsigned int wait_loc; // for packets waiting to be injected into the network
    unsigned int wait_type;

    WARPED_REGISTER_SERIALIZABLE_MEMBERS(cereal::base_class<warped::Event>(this),
                                            receiver_name_, event_ts_)
};


class Node : public warped::LogicalProcess {
public:
    Node( const string& name,
            unsigned int grid_dimension,
            unsigned int grid_size,
            unsigned int grid_order,
            unsigned int dimension_index,
            unsigned int index )

        :   LogicalProcess(name),
            state_(),
            grid_dimension_(grid_dimension),
            grid_size_(grid_size),
            grid_order_(grid_order),
            dimension_index_(dimension_index),
            index_(index) {
    }

    virtual warped::LPState& getState() { return state_; }

    virtual vector<shared_ptr<warped::Event> > initializeLP() override;

    virtual vector<shared_ptr<warped::Event> > receiveEvent(const warped::Event&);

    NodeState state_;

    unsigned int grid_dimension_;
    unsigned int grid_size_;
    unsigned int grid_order_;
    unsigned int dimension_index_;
    unsigned int index_;

protected:
    unsigned int neighbor(unsigned int index, unsigned int destination);
};

#endif
