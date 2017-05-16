//Implementation of a Forest Fire Simulation
#ifndef FOREST_HPP_DEFINED
#define FOREST_HPP_DEFINED

#include <string>
#include <vector>
#include <memory>
#include <random>

#include "warped.hpp"


WARPED_DEFINE_LP_STATE_STRUCT(ForestState) {
    unsigned int burning_status_;  //Range from 0 to 5 respresenting the amount of vegitation burned in an LP
};

enum forest_event_t {
    IGNITION,       //Another LP spreads and starts a fire in this LP
    SPREADING       //The fire in this current LP is spreading to another LP
};

enum direction_t {  //The direction of the spread of fire from the currently burning LP

    NORTH,
    NORTH_WEST,
    NORTH_EAST,
    SOUTH,
    SOUTH_WEST,
    SOUTH_EAST,
    EAST,
    WEST,
};

class ForestEvent : public warped::Event {  //Class Definition of a forest fire event
    public:
        ForestEvent() = default;
        ForestEvent(const std::string& receiver_name, const forest_event_t type,
                                                        const unsigned int timestamp)
            : receiver_name_(receiver_name), type_(type), ts_(timestamp) {} //Initializing Class Components
        
        const std::string& receiverName() const { return receiver_name_; }
        unsigned int timestamp() const { return ts_; }

        std::string receiver_name_;
        forest_event_t type_;
        unsigned int ts_;

        WARPED_REGISTER_SERIALIZABLE_MEMBERS(cereal::base_class<warped::Event>(this), receiver_name_, type_, ts_)

};

class Forest : public warped::LogicalProcess{
    public:
        Forest( const std::string& name,
                const unsigned int num_forest_pixel_x,    
                const unsigned int num_forest_pixel_y,    
                const unsigned int vegetation_type, 
                const unsigned int elevation,   
                const unsigned int moisture,    
                const unsigned int ignition_mean,   
                const unsigned int spread_mean,
                const unsigned int index)
                
        :   LogicalProcess(name),
            state_(),
            rng_(new std::default_random_engine(index)),
            num_forest_pixel_x_(num_forest_pixel_x),    
            num_forest_pixel_y_(num_forest_pixel_y),    
            vegetation_type_(vegetation_type), 
            elevation_(elevation),   
            vegetation_moisture_(vegetation_moisture),    
            ignition_mean_(ignition_mean),   
            spread_mean_(spread_mean),
            index_(index){
                state_.burning_status_ = 0; //Should we create this state as an enum instead of an unsigned int????????????
            }
            
            virtual std::vector<std::shared_ptr<warped::Event> > initializeLP() override;
            virtual std::vector<std::shared_ptr<warped::Event> > receiveEvent(const warped::Event&);
            virtual warped::LPState& getState() { return this->state_; }
            
            ForestState state_;
            
            static inline std::string lp_name(const unsigned int);
            
            
            
            
    protected:
        std::shared_ptr<std::default_random_engine> rng_;
        const unsigned int num_forest_pixel_x_;   //the x coordinate of the LP on the picture
        const unsigned int num_forest_pixel_y_;    //the y coordinate of the LP on the picture
        const unsigned int vegetation_type_; //the type of vegetation growing in the LP
        const unsigned int elevation_;   //the elevation of the LP
        const unsigned int vegetation_moisture_;    //the amount of moisture in the LP
        const unsigned int ignition_mean_;   
        const unsigned int spread_mean_;
        const unsigned int index_; //The identifier used by the model to distinguish between LPs
        
        std::string compute_move(direction_t direction);

}

#endif