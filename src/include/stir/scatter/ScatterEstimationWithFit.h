#ifndef __scatter_ScatterEstimationWithFit_H__
#define __scatter_ScatterEstimationWithFit_H__

#include "stir/shared_ptr.h"
#include "stir/scatter/ScatterEstimationByBin.h"

START_NAMESPACE_STIR

class ScatterEstimationWithFit : public ScatterEstimationByBin
{
private:
    typedef ScatterEstimationByBin base_type;

public:
    ScatterEstimationWithFit();

    virtual Succeeded process_data();



protected:
    void set_defaults();
    void initialise_keymap();
    bool post_processing();

    //!
    //! \brief pre_processing
    //! \return
    //! \details This function will check if the proper images are supplied.
    //! Will make the zoom, normalization calculations etc.
    bool pre_processing();

    // for density_image_for_scatter_points_filename;
    // todo --- reconsider them later.
    //!
    //! \brief transverse_scatter_point_sampling_in_mm
    //! \details Transverse scatter point sampling distance
    float transverse_scatter_point_sampling_in_mm;
    //!
    //! \brief axial_scatter_point_sampling_in_mm
    //! \details Axial scatter point sampling distance
    float axial_scatter_point_sampling_in_mm;

    // for template_proj_data_filename
    // - If not initialized otherwise
    //!
    //! \brief num_detectors_per_ring_for_scatter_simulation
    //! \details The number of detectors per ring in the subsampled scanner template
    int num_detectors_per_ring_for_scatter_simulation;

    //!
    //! \brief num_detectors_rings_for_scatter_simulation
    //! \details The number of rings in the subsampled scanner template
    int num_detectors_rings_for_scatter_simulation;

};
END_NAMESPACE_STIR

#include "stir/scatter/ScatterEstimationWithFit.inl"

#endif
