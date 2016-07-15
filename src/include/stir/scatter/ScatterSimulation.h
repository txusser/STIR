#ifndef __stir_scatter_ScatterSimulation_H__
#define __stir_scatter_ScatterSimulation_H__

/*!
  \file
  \ingroup scatter
  \brief Definition of class stir::ScatterSimulation

  \author Nikos Efthimiou
  \author Kris Thielemans
*/

#include "stir/ParsingObject.h"
#include "stir/Succeeded.h"
#include "stir/RegisteredObject.h"

#include "stir/VoxelsOnCartesianGrid.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"

START_NAMESPACE_STIR

/*!
  \ingroup scatter
  \brief Simuate the scatter probability using a model-based approach

  N.E. : This class is roughly the base class of what used to be the ScatterEstimationByBin.
  Because there are different approaches on the actual simulation process, this base class will
  be in charge of hold projection data and subsample the attenuation image, while more core function
  will deligate to classes like SignelScatterSimulation.

  This class computes the single Compton scatter estimate for PET data using an analytical
  approximation of an integral. It takes as input an emission image and an attenuation image.
  This is effectively an implementation of the simulation part of the algorithms of
  Watson and Ollinger.

  One non-standard feature is that you can specify a different attenuation image to find the
  scatter points and one to compute the integrals over the attenuation image. The idea is that
  maybe you want to compute the integrals on a finer grid than you sample the attenuation image.
  This is probably not very useful though.

  \todo Currently this can only be run by initialising it via parsing of a file. We need
  to add a lot of set/get members.

  \todo This class currently uses a simple Gaussian model for the energy resolution. This model
  and its parameters (\a reference_energy and \a energy_resolution) really should be moved the
  the Scanner class. Also the \a lower_energy_threshold and \a upper_energy_threshold should
  be read from the emission data, as opposed to setting them here.

  \todo detector coordinates are derived from ProjDataInfo, but areas and orientations are
  determined by using a cylindrical scanner.

  \todo This class should be split into a generic class and one specific to PET single scatter.

  \par References
  This implementation is described in the following
  <ol>
  <li> C. Tsoumpas, P. Aguiar, K. S. Nikita, D. Ros, K. Thielemans,
  <i>Evaluation of the Single Scatter Simulation Algorithm Implemented in the STIR Library</i>,
  Proc. of IEEE Medical Imaging Conf. 2004, vol. 6 pp. 3361 - 3365.
  </li>
  </ol>
  Please refer to the above if you use this implementation.

  See also these papers for extra evaluation

  <ol>
  <li>N. Dikaios , T. J. Spinks , K. Nikita , K. Thielemans,
     <i>Evaluation of the Scatter Simulation Algorithm Implemented in STIR,</i>
     proc. 5th ESBME, Patras, Greece.
  </li>
  <li>P. Aguiar, Ch. Tsoumpas, C. Crespo, J. Pavia, C. Falcon, A. Cot, K. Thielemans and D. Ros,
     <i>Assessment of scattered photons in the quantification of the small animal PET studies,</i>
     Eur J Nucl Med Mol I 33:S315-S315 Sep 2006, Proc. EANM 2006, Athens, Greece.
  </li>
  </ol>
*/

class ScatterSimulation : public RegisteredObject<ScatterSimulation>,
        public ParsingObject
{

public:

    //!
    //! \brief ScatterSimulation
    //! \details Default constructor
    ScatterSimulation();

    virtual ~ScatterSimulation();

    virtual Succeeded
    process_data();

    //! gives method information
    virtual std::string method_info() const = 0;

    //! prompts the user to enter parameter values manually
    virtual void ask_parameters();

protected:

    void initialise(const std::string& parameter_filename);

    virtual void set_defaults();

    virtual void initialise_keymap();

    //! used to check acceptable parameter ranges, etc...
    virtual bool post_processing();

    //!
    //! \brief sub_atten_image_sptr
    //! \details This is the image with tha anatomical information
    //! \warning Previously density_image_for_scatter_points_sptr
    shared_ptr<VoxelsOnCartesianGrid<float> > sub_atten_image_sptr;

    //!
    //! \brief sub_proj_data_info_ptr
    //!
    ProjDataInfoCylindricalNoArcCorr * sub_proj_data_info_ptr;

    //! computes scatter for one viewgram
    /*! \return total scatter estimated for this viewgram */
    virtual double
      process_data_for_view_segment_num(const ViewSegmentNumbers& vs_num);

private:

    /**
      *
      * \name Variables related to the subsampled attenuation image.
      * @{
    */

    //!
    //! \brief recompute_sub_atten_image
    bool recompute_sub_atten_image;

    //!
    //! \brief sub_atten_image_filename
    //! \details Input or output file name of the subsampled
    //! attenuation image, depends on the reconmpute_sub_atten_image
    std::string sub_atten_image_filename;

    //!
    //! \brief sub_vox_xy
    //! \details The subsampling of the attenuation image is done,
    //! in the arbitary zoom factors. This correspond to the zoom in
    //! the XY plane.
    float sub_vox_xy;

    //!
    //! \brief sub_vox_z
    //! \details The subsampling of the attenuation image is done,
    //! in the arbitary zoom factors. This correspond to the zoom in
    //! the Z axis.
    float sub_vox_z;

    /** }@*/

    /**
     *
     * \name Variables related to the subsampled ProjData.
     * @{
     */

    //!
    //! \brief recompute_sub_projdata
    //! \details
    bool recompute_sub_projdata;

    //!
    //! \brief sub_template_proj_data_filename
    //!
    std::string sub_proj_data_filename;

    //!
    //! \brief sub_num_dets_per_ring
    //! \details Number of detectors per ring for the subsampled
    //! scanner
    int sub_num_dets_per_ring;

    //!
    //! \brief sub_num_rings
    //! \details Number of rings for the subsampled scanner.
    int sub_num_rings;

    /** }@*/


    //!
    //! \brief subsample_image
    //! \param _this_image_sptr The image which will be subsampled
    //! \param _new_image_sptr The new subsampled image
    //! \param _sub_vox_xy The zoom factor on the xy plane
    //! \param _sub_vox_z The zoom factor on the z axis
    //! \param output_filename If set, the file to store the subsampled image
    //! \details Replaces the zoom_att_image.sh.
    void subsample_image(shared_ptr<VoxelsOnCartesianGrid<float> > & _this_image_sptr,
                         shared_ptr<VoxelsOnCartesianGrid<float> > & _new_image_sptr,
                         float& _sub_vox_xy, float& _sub_vox_z,
                         std::string output_filename = "");

    //!
    //! \brief subsample_projdata_info
    //! \param _original_projdata_info
    //! \param _new_projdata_info
    //! \param output_filename
    //! \details A function to create a subsampled scanner and projdata info
    //! for a new number of detectors and rings.
    void subsample_projdata_info(ProjDataInfoCylindricalNoArcCorr* _original_projdata_info,
                            int new_num_dets_per_ring, int new_num_rings,
                            std::string output_filename = "");


};

END_NAMESPACE_STIR

#endif


