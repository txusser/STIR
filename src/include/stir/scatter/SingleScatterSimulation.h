/*!
  \file
  \ingroup scatter
  \brief Definition of class stir::SingleScatterSimulation

  \author Nikos Efthimiou
  \author Kris Thielemans
*/

#ifndef __stir_scatter_SingleScatterSimulation_H__
#define __stir_scatter_SingleScatterSimulation_H__

#include "stir/Succeeded.h"
#include "stir/scatter/ScatterSimulation.h"
#include "stir/shared_ptr.h"
#include "stir/RegisteredParsingObject.h"


START_NAMESPACE_STIR

class SingleScatterSimulation : public
        RegisteredParsingObject<
        SingleScatterSimulation,
        ScatterSimulation,
        ScatterSimulation >
{
private:
    typedef RegisteredParsingObject<
    SingleScatterSimulation,
    ScatterSimulation,
    ScatterSimulation > base_type;
public:

    //! Name which will be used when parsing a OSMAPOSLReconstruction object
    static const char * const registered_name;

    //!
    //! \brief ScatterSimulation
    //! \details Default constructor
    SingleScatterSimulation();

    //!
    //! \brief ScatterSimulation
    //! \param parameter_filename
    //! \details Constructor with initialisation from parameter file
    explicit
    SingleScatterSimulation(const std::string& parameter_filename);

    virtual ~SingleScatterSimulation();

    //    virtual Succeeded
    //    process_data();

    //! gives method information
    virtual std::string method_info() const;

    //! prompts the user to enter parameter values manually
    virtual void ask_parameters();

protected:

    void initialise(const std::string& parameter_filename);

    virtual void set_defaults();
    virtual void initialise_keymap();

    //! used to check acceptable parameter ranges, etc...
    virtual bool post_processing();


};

END_NAMESPACE_STIR

#endif
