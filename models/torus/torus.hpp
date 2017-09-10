#ifndef TORUS_HPP
#define TORUS_HPP

#include <string>
#include <memory>
#include <random>

#include "warped.hpp"

WARPED_DEFINE_LP_STATE_STRUCT(NodeState) {
    bool available_ = true;
};

class MessageEvent : public warped::Event {
public:
    MessageEvent() = default;

    MessageEvent( const std::string   receiver_name,
                  unsigned int        event_ts  )

        :   receiver_name_(receiver_name),
            event_ts_(event_ts) {}

    const std::string& receiverName() const { return receiver_name_; }

    unsigned int timestamp() const { return event_ts_; }

    unsigned int size() const {
        return receiver_name_.length() + sizeof(event_ts_);
    }

    std::string     receiver_name_;
    std::string     destination_name;
    unsigned int  destination;
    unsigned int    hop_count_;
    unsigned int    event_ts_;

    WARPED_REGISTER_SERIALIZABLE_MEMBERS(cereal::base_class<warped::Event>(this),
                                            receiver_name_, event_ts_)
};


class Node : public warped::LogicalProcess {
public:
    Node( const std::string& name,
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

    virtual std::vector<std::shared_ptr<warped::Event> > initializeLP() override;

    virtual std::vector<std::shared_ptr<warped::Event> > receiveEvent(const warped::Event&);

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
