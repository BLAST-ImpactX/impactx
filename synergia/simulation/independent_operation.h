#ifndef INDEPENDENT_OPERATION_H_
#define INDEPENDENT_OPERATION_H_

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/utils/logger.h"

class Independent_operation {
private:
  std::string type;

private:
  virtual void
  print_impl(Logger& logger) const
  {}
  virtual void apply_impl(Bunch& bunch, Logger& logger) const = 0;

public:
  Independent_operation(std::string const& type) : type(type) {}
  virtual ~Independent_operation() = default;

  void
  apply(Bunch& bunch, Logger& logger) const
  {
    apply_impl(bunch, logger);
  }

  std::string const&
  get_type() const
  {
    return type;
  }

  void
  print(Logger& logger) const
  {
    logger(LoggerV::INFO_OPN) << "\ttype = " << type << ", ";
    print_impl(logger);
  }
};

class LibFF_operation : public Independent_operation {
private:
  std::vector<Lattice_element_slice> slices;

private:
  void
  print_impl(Logger& logger) const override
  {}
  void apply_impl(Bunch& bunch, Logger& logger) const override;

public:
  LibFF_operation(std::vector<Lattice_element_slice> const& slices)
    : Independent_operation("LibFF"), slices(slices)
  {}
};

#endif /* INDEPENDENT_OPERATION_H_ */
