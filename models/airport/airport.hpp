// An implementation of Fujimoto's airport model
// Ported from the ROSS airport model (https://github.com/carothersc/ROSS/blob/master/ross/models/airport)
// Author: Eric Carver (carverer@mail.uc.edu)

#ifndef AIRPORT_HPP_DEFINED
#define AIRPORT_HPP_DEFINED

#include <string>
#include <vector>
#include <memory>

#include "warped.hpp"

#include "MLCG.h"

WARPED_DEFINE_OBJECT_STATE_STRUCT(AirportState) {
  unsigned int landings_;
  unsigned int departures_;
  unsigned int planes_flying_;
  unsigned int planes_grounded_;
};

enum airport_event_t {
  ARRIVAL,
  DEPARTURE,
  LANDING
};

enum direction_t {
  NORTH,
  EAST,
  SOUTH,
  WEST
};

class AirportEvent : public warped::Event {
public:
  AirportEvent() = default;
  AirportEvent(const std::string& receiver_name, const airport_event_t type, 
                const unsigned int timestamp)
    : receiver_name_(receiver_name), type_(type), ts_(timestamp) {}

  const std::string& receiverName() const { return receiver_name_; }
  unsigned int timestamp() const { return ts_; }

  std::string receiver_name_;
  airport_event_t type_;
  unsigned int ts_;

  WARPED_REGISTER_SERIALIZABLE_MEMBERS(cereal::base_class<warped::Event>(this), receiver_name_, type_, ts_)
};

class Airport : public warped::SimulationObject {
public:
  Airport(const std::string& name, const unsigned int num_airports, const unsigned int num_planes,
          const unsigned int land_mean, const unsigned int depart_mean, const unsigned int index)
    : SimulationObject(name), state_(), rng_(new MLCG), num_airports_(num_airports), num_planes_(num_planes),
      land_mean_(land_mean), depart_mean_(depart_mean), index_(index) {}

  virtual std::vector<std::shared_ptr<warped::Event> > createInitialEvents();
  virtual std::vector<std::shared_ptr<warped::Event> > receiveEvent(const warped::Event&);
  virtual warped::ObjectState& getState() { return this->state_; }

  AirportState state_;

  static inline std::string object_name(const unsigned int);

protected:
  std::shared_ptr<MLCG> rng_;
  std::default_random_engine rng_engine_;
  const unsigned int num_airports_;
  const unsigned int num_planes_;
  const unsigned int land_mean_;
  const unsigned int depart_mean_;
  const unsigned int index_;
};

#endif
