#ifndef FIXED_T_Z_CONVERTER_H_

// #include "synergia/utils/cereal.h"
class Bunch;

/// Fixed_t_z_converter is a virtual base class for converters
/// between the different system of coordinates

/// the collective effects addressed with the split-operator method (Yoshida)
/// should always to be done in the canonical coordinates (x, px, y, py, t, -E)

/// space-charge and wake field calculation always assume a rigid beam
/// approximation i.e. all the particles move with the same velocity v=beta c,
/// rho(x,y,z,t)=rho(x,y,z-vt) for these cases, a transformation of coordinates
/// which assume all particles move longitudinally with the same beta c (i.e.
/// beta_x_i=beta_y_i=0 and beta_z_i = beta) is appropriate converters to and
/// from bunch frame will be removed once the space_charge_hockney is revisited

class Fixed_t_z_converter {
  public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
    /// state in the accelerator frame.
    virtual void from_z_lab_to_t_lab(Bunch& bunch) = 0;
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z
    /// state in the accelerator frame.
    virtual void from_t_lab_to_z_lab(Bunch& bunch) = 0;
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
    /// state in the beam frame. Should never be used!!!!!
    virtual void from_z_lab_to_t_bunch(Bunch& bunch) = 0;
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in
    /// the accelerator frame. Should never be used!!!!!
    virtual void from_t_bunch_to_z_lab(Bunch& bunch) = 0;
    /// Should never be used!!!!!
    virtual void from_t_lab_to_t_bunch(Bunch& bunch) = 0;
    /// Should never be used!!!!!
    virtual void from_t_bunch_to_t_lab(Bunch& bunch) = 0;

    //
    //     virtual void
    //     fixed_z_lab_to_z_bunch(Bunch &bunch) = 0;
    //

    /// Convert from the fixed-t state to the fixed-z state.
    // virtual void
    // fixed_t_to_fixed_z(Bunch &bunch) = 0;

    /// Convert from the fixed-z state to the fixed-t state.
    //  virtual void
    //  fixed_z_to_fixed_t(Bunch &bunch) = 0;

#if 0
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
#endif

    virtual ~Fixed_t_z_converter(){};
};
#if 0
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_converter)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(Fixed_t_z_converter);
#endif

/// Fixed_t_z_zeroth implements a fixed-t-fixed-z converter using
/// the simplest approximation: longitudinal coordinates are transformed,
/// but transverse coordinates are unaffected.
class Fixed_t_z_zeroth : public Fixed_t_z_converter {
  public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
    /// state in the accelerator frame.
    void from_z_lab_to_t_lab(Bunch& bunch);
    /// Convert from the fixed-t state in the accelerator frameto the fixed-z
    /// state in the accelerator frame.
    void from_t_lab_to_z_lab(Bunch& bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
    /// state in the beam frame. Should never be used!!!!!
    void from_z_lab_to_t_bunch(Bunch& bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in
    /// the accelerator frame. Should never be used!!!!!
    void from_t_bunch_to_z_lab(Bunch& bunch);
    /// Should never be used!!!!!
    void from_t_lab_to_t_bunch(Bunch& bunch);
    /// Should never be used!!!!!
    void from_t_bunch_to_t_lab(Bunch& bunch);

#if 0
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
#endif
};
#if 0
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_zeroth)
#endif

/// Fixed_t_z_zeroth implements a fixed-t-fixed-z converter using
/// the ballistic approximation: longitudinal coordinates are transformed,
/// then transverse coordinates are transformed using the ballistic
/// approximation, i.e., as though the particles were traveling independently
/// in free space.
/// UNIMPLEMENTED
/*class Fixed_t_z_ballistic : public Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
state in the accelerator frame. void from_z_lab_to_t_lab(Bunch &bunch)
    {
    }
    ;
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z
state in the accelerator frame. void from_t_lab_to_z_lab(Bunch &bunch)
    {
    }
    ;
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
state in the beam frame. void from_z_lab_to_t_bunch(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in
the accelerator frame. void from_t_bunch_to_z_lab(Bunch &bunch); void
    from_t_lab_to_t_bunch(Bunch &bunch){};
    void
    from_t_bunch_to_t_lab(Bunch &bunch){};

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_ballistic)

class Fixed_t_z_alex : public Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
state in the accelerator frame. void from_z_lab_to_t_lab(Bunch &bunch);
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z
state in the accelerator frame. void from_t_lab_to_z_lab(Bunch &bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
state in the beam frame. void from_z_lab_to_t_bunch(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in
the accelerator frame. void from_t_bunch_to_z_lab(Bunch &bunch); void
     from_t_lab_to_t_bunch(Bunch &bunch);

     void
     from_t_bunch_to_t_lab(Bunch &bunch);

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_alex)

/// transformation as in the old synergia....
class Fixed_t_z_synergia20 : public Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
state in the accelerator frame. void from_z_lab_to_t_lab(Bunch &bunch);
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z
state in the accelerator frame. void from_t_lab_to_z_lab(Bunch &bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t
state in the beam frame. void from_z_lab_to_t_bunch(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in
the accelerator frame. void from_t_bunch_to_z_lab(Bunch &bunch); void
    from_t_lab_to_t_bunch(Bunch &bunch);
    void
    from_t_bunch_to_t_lab(Bunch &bunch);

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_synergia20)*/

#define FIXED_T_Z_CONVERTER_H_

#endif /* FIXED_T_Z_CONVERTER_H_ */
