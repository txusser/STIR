#ifndef __stir_scatter_ScatterEstimationByBin_H__
#define __stir_scatter_ScatterEstimationByBin_H__

/*
    Copyright (C) 2004 - 2009 Hammersmith Imanet Ltd
    Copyright (C) 2013 University College London
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    See STIR/LICENSE.txt for details
*/
/*!
  \file
  \ingroup scatter
  \brief Definition of class stir::ScatterEstimationByBin.
  
  \author Charalampos Tsoumpas
  \author Nikolaos Dikaios
  \author Kris Thielemans
*/

#include "stir/shared_ptr.h"
#include "stir/DiscretisedDensity.h"
#include "stir/ProjData.h"
#include "stir/ParsingObject.h"
#include "stir/numerics/BSplines.h"
#include <vector>
#include "stir/CartesianCoordinate3D.h"
#include "stir/recon_buildblock/Reconstruction.h"

#include "stir/VoxelsOnCartesianGrid.h"
#include "stir/IndexRange3D.h"

#include "stir/scatter/ScatterSimulation.h"

START_NAMESPACE_STIR

class Succeeded;
class ProjDataInfoCylindricalNoArcCorr;
class ViewSegmentNumbers;
class BinNormalisation;

/*!
  \ingroup scatter
  \brief Estimate the scatter probability using a model-based approach

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
class ScatterEstimationByBin : public ParsingObject
{
 public:
  //! upsample coarse scatter estimate and fit it to tails of the emission data
  /*! Current procedure:
    1. interpolate segment 0 of \a scatter_proj_data to size of segment 0 of \a emission_proj_data
    2. inverseSSRB to create oblique segments
    3. find scale factors with get_scale_factors_per_sinogram()
    4. apply thresholds
    5. filter scale-factors in axial direction (independently for every segment)
    6. apply scale factors using scale_sinograms()
  */
 static void
   upsample_and_fit_scatter_estimate(ProjData& scaled_scatter_proj_data,
                     const  ProjData& emission_proj_data,
                     const ProjData& scatter_proj_data,
                                     const BinNormalisation& scatter_normalisation,
                     const ProjData& weights_proj_data,
                     const float min_scale_factor,
                     const float max_scale_factor,
                     const unsigned half_filter_width,
                     BSpline::BSplineType spline_type,
                     const bool remove_interleaving = true);


  //! Default constructor (calls set_defaults())
  ScatterEstimationByBin();

  virtual Succeeded process_data();


  // N.E: New main process_data funtion,
  // the '_' will be removed later.
  virtual Succeeded _process_data();

  virtual Succeeded
  _iterate(int,
           shared_ptr<ExamData>&,
           shared_ptr<ExamData>&,
           shared_ptr<ExamData>&,
           shared_ptr<VoxelsOnCartesianGrid<float> >&);

  /**
   *
   *  \name functions to (re)set images or projection data
   *  These functions also invalidate cached activity integrals such that the cache will be recomputed.
   *
   *  The functions that read a file call error() if the reading failed.
   *
   * @{
   */

  //!
  //! \brief set_activity_image_sptr
  //!
  inline Succeeded set_activity_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >&);

  //!
  //! \brief set_activity_image
  //! \param filename
  //!
  inline void set_activity_image(const std::string& filename);

  inline Succeeded set_atten_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >&);

  //!
  //! \brief set_image_from_file
  //! \param filename
  //! \details This function loads image files from the disk
  inline shared_ptr<DiscretisedDensity<3,float> >
    get_image_from_file(const std::string& filename);

  //!
  //! \brief set_template_proj_data_info_from_file
  //! \param filename
  //! \details Sets the template of the ProjDataInfo from a file
  inline void set_template_proj_data_info_from_file(const std::string& filename);

  //!
  //! \brief set_sub_proj_data_info_from_file
  //! \param filename
  //!
  inline void
  set_sub_proj_data_info_from_file(const std::string& filename);

  //!
  //! \brief set_atten_coeffs_from_file
  //! \param filename
  //! \details Load attenuation coefficients from disk, instead of
  //! calculating them on the fly.
  //! \warning Not tested!
  inline void set_atten_coeffs_from_file(const std::string& filename);

  inline void set_template_proj_data_info_sptr(const shared_ptr<ProjDataInfo>&);

  //! set the image that determines where the scatter points are
  /*! Also calls sample_scatter_points()
   \warning Uses attenuation_threshold member variable
  */
  inline Succeeded set_sub_atten_image_sptr(const shared_ptr<DiscretisedDensity<3,float> >&);


  //void set_output_proj_data_sptr(const shared_ptr<ProjData>& new_sptr);
  //! create output projection data of same size as template_proj_data_info
  /*! \warning use set_template_proj_data_info() first.

   Currently always uses Interfile output.
   \warning If the specified file already exists it will be erased.
  */
  inline void set_proj_data_from_file(const std::string& filename,
                                      shared_ptr<ProjData>& _this_projdata);

  /** @}*/

  // TODO write_log can't be const because parameter_info isn't const
  virtual void
    write_log(const double simulation_time,
          const float total_scatter);


  /**
   *
   * \name functions to set parameters
   * @{
   */

    inline void set_attenuation_threshold(float _val);

    inline void set_random_point(bool _val);

    inline void set_cache_enabled(bool _val);

  /** @}*/


 protected:
  void set_defaults();
  void initialise_keymap();
  bool post_processing();

  //! threshold below which a voxel in the attenuation image will not be considered as a candidate scatter point
  float attenuation_threshold;

  //! boolean to see if we need to move the scatter point randomly within in its voxel
  /*! This was first recommended by Watson. It is recommended to leave this on, as otherwise
     discretisation errors are more obvious.

     Note that the random generator is seeded via date/time, so re-running the scatter
     simulation will give a slightly different result if this boolean is on.
  */
  bool random;
  //! boolean to see if we need to cache the integrals
  /*! By default, we cache the integrals over the emission and attenuation image. If you run out
      of memory, you can switch this off, but performance will suffer dramatically.
  */
  bool use_cache;

  /**
    *
    * \name Variables retated to the scanner geometry
    * @{
    */

  //!
  //! \brief template_proj_data_filename
  //! \details The file name for the Scanner template
  std::string template_proj_data_filename;

  //!
  //! \brief proj_data_info_ptr
  //! \details The projection data info of the scanner template
  ProjDataInfoCylindricalNoArcCorr * proj_data_info_ptr;


  //!
  //! \brief proj_data_info_2d_ptr
  //! \details The projection data info for the 2D data - after SSRB
  ProjDataInfoCylindricalNoArcCorr * proj_data_info_2d_ptr;

  //!
  //! \brief proj_data_info_2d_sptr
  //! \details shared pointer to projection data info for the 2D data - after SSRB
  shared_ptr < ProjDataInfo > proj_data_info_2d_sptr;

  //!
  //! \brief template_exam_info_sptr
  //! \details Exam info extracted from the scanner template
  shared_ptr<ExamInfo> template_exam_info_sptr;

  /** @}*/


  /**
   *
   * \name Variables related to the initial activity image and the measured emission data.
   * @{
   */

  //!
  //! \brief recompute_initial_estimate
  //! \details If set then the initial activity estimate will be recomputed
  //! and stored if a name is provided.
  bool recompute_initial_activity_image;

  //!
  //! \brief initial_activity_image_filename
  //! \details Filename of the initial activity image.
  std::string initial_activity_image_filename;

  //!
  //! \brief reconstruction_method_sptr
  //! \details The reconsturction which is going to be used for the scatter simulation
  //! and the intial activity image (if recompute set).
  shared_ptr < Reconstruction < DiscretisedDensity < 3, float > > >
          reconstruction_template_sptr;

  //!
  //! \brief reconstruction_template_par_filename
  //! \details The filename for the parameters file of the reconstruction method.
  //! \warning Refere to the samples for a proper example.
  std::string reconstruction_template_par_filename;

  //!
  //! \brief activity_image_sptr
  //! \details Initially with is the reconstructed activity image, but during the scatter
  //! estimation it with actually hold the iterative estimates.
  //! Therefore the nane might change later.
  shared_ptr<VoxelsOnCartesianGrid<float> > activity_image_sptr;

  //!
  //! \brief input_projdata_filename
  //! \details Filename of the measured emission data.
  std::string input_projdata_filename;

  //!
  //! \brief input_data
  //! \details This memnbers holds the measured emission data.
  shared_ptr<ProjData> input_projdata_sptr;

  shared_ptr<ProjData> input_projdata_ssrb_sptr;
  /** }@*/



  /**
   *
   * \name Variables related to the attenuation image and coefficients
   * @{
   */

  //!
  //! \brief atten_image_filename
  //! \details This is the image file name with the anatomic information.
  //! \warning Previously named density_image.
  std::string atten_image_filename;

  //!
  //! \brief recompute_atten_coeff
  //! \details If set to 1 the attenuation coefficients are going to
  //! be recalculated.
  bool recompute_atten_coeff;

  //!
  //! \brief atten_coeff_filename
  //! \details The file name for the attenuation coefficients.
  //! If they are to be recalculated they will be stored here, if set.
  std::string atten_coeff_filename;

  //!
  //! \brief atten_image_sptr
  //!
  shared_ptr<VoxelsOnCartesianGrid<float> > atten_image_sptr;

  //!
  //! \brief output_proj_data_sptr
  //! \details The attenuation coefficients
  shared_ptr<ProjData> atten_coeffs_sptr;

  /** @}*/

  /**
    * \name Varianbles realted to the background proj data
    * @{
    *
    */

  std::string back_projdata_filename;

  shared_ptr<ProjData> back_projdata_sptr;

  /** @}*/

  /**
   *
   *  \name Variables related to normalization factors
   *  @{
   */

  //!
  //! \brief normalization_coeffs_filename
  //! \details File with normalization factors
  std::string normalization_coeffs_filename;

  //!
  //! \brief normalization_factors_sptr
  //! \details Sinogram with the normalization factors.
  shared_ptr<ProjData> normalization_factors_sptr;

  /** @}*/

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
  //! \brief sub_atten_image_sptr
  //! \details This is the image with tha anatomical information
  //! \warning Previously density_image_for_scatter_points_sptr
  shared_ptr<VoxelsOnCartesianGrid<float> > sub_atten_image_sptr;

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


  std::string mask_image_filename;

  std::string mask_parameters_filename;

  bool recompute_mask_image;

  shared_ptr < DiscretisedDensity < 3, float >  > mask_image_sptr;

  float mask_max_threshold;
  float mask_add_scalar;
  float mask_min_threshold;
  float mask_times_scalar;

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
  //! \brief sub_proj_data_info_ptr
  //!
  ProjDataInfoCylindricalNoArcCorr * sub_proj_data_info_ptr;

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

  std::string output_proj_data_filename;
  shared_ptr<ProjData> output_proj_data_sptr;


  /*************** functions that do the work **********/

  enum image_type{act_image_type, att_image_type};
  struct ScatterPoint
  {
    CartesianCoordinate3D<float> coord;
    float mu_value;
  };

  std::vector< ScatterPoint> scatt_points_vector;
  float scatter_volume;

  //! find scatter points
  /*! This function sets scatt_points_vector and scatter_volume. It will also
      remove any cached integrals as they would be incorrect otherwise.
  */
  void
    sample_scatter_points();


  /************************************************************************/

  /** \name detection related functions
   *
   * @{
   */

  //!
  //! \brief detection_efficiency
  //! \param energy
  //! \return
  //! \details energy-dependent detection efficiency (Gaussian model)
  float detection_efficiency(const float energy) const;


  //! maximum angle to consider above which detection after Compton scatter is considered too small
  static
    float
    max_cos_angle(const float low, const float approx, const float resolution_at_511keV);

  //! mimumum energy to consider above which detection after Compton scatter is considered too small
  static
    float
    energy_lower_limit(const float low, const float approx, const float resolution_at_511keV);

  virtual
    void
    find_detectors(unsigned& det_num_A, unsigned& det_num_B, const Bin& bin) const;

  unsigned
    find_in_detection_points_vector(const CartesianCoordinate3D<float>& coord) const;
  // private:


  CartesianCoordinate3D<float>  shift_detector_coordinates_to_origin;

  //! average detection efficiency of unscattered counts
  double
    detection_efficiency_no_scatter(const unsigned det_num_A,
                    const unsigned det_num_B) const;

  // next needs to be mutable because find_in_detection_points_vector is const
  mutable std::vector<CartesianCoordinate3D<float> > detection_points_vector;
 private:
  int total_detectors;

  //@}
 protected:
  //! computes scatter for one viewgram
  /*! \return total scatter estimated for this viewgram */
  virtual double
    process_data_for_view_segment_num(const ViewSegmentNumbers& vs_num);

  /************************************************************************/

  //! \name integrating functions
  //@{
  static
    float
    integral_between_2_points(const DiscretisedDensity<3,float>& density,
                  const CartesianCoordinate3D<float>& point1,
                  const CartesianCoordinate3D<float>& point2);

  float
    exp_integral_over_attenuation_image_between_scattpoint_det (const CartesianCoordinate3D<float>& scatter_point,
                                const CartesianCoordinate3D<float>& detector_coord);



  float
    integral_over_activity_image_between_scattpoint_det (const CartesianCoordinate3D<float>& scatter_point,
                             const CartesianCoordinate3D<float>& detector_coord);


    
  float
    cached_integral_over_activity_image_between_scattpoint_det(const unsigned scatter_point_num,
                                   const unsigned det_num);

  float
    cached_exp_integral_over_attenuation_image_between_scattpoint_det(const unsigned scatter_point_num,
                                      const unsigned det_num);
  //@}

  /************************************************************************/

  float
    single_scatter_estimate_for_one_scatter_point(const std::size_t scatter_point_num,
                       const unsigned det_num_A,
                       const unsigned det_num_B);


  void
    single_scatter_estimate(double& scatter_ratio_singles,
                const unsigned det_num_A,
                const unsigned det_num_B);


  virtual double
    scatter_estimate(const unsigned det_num_A,
             const unsigned det_num_B);



 public:
  /************************************************************************/

  //! \name Compton scatter cross sections
  //@{
  static
    inline float
    dif_Compton_cross_section(const float cos_theta, float energy);


  static
    inline float
    total_Compton_cross_section(float energy);

  static
    inline float
    photon_energy_after_Compton_scatter(const float cos_theta, const float energy);

  static
    inline float
    photon_energy_after_Compton_scatter_511keV(const float cos_theta);

  static
    inline float
    total_Compton_cross_section_relative_to_511keV(const float energy);
  //@}
  /************************************************************************/
 protected:

  float
    compute_emis_to_det_points_solid_angle_factor(const CartesianCoordinate3D<float>& emis_point,
                          const CartesianCoordinate3D<float>& detector_coord) ;

  //!
  //! \brief calculate_attenuation_coefficients
  //! \param proj_data_info_ptr
  //! \param atten_image_sptr
  //! \return
  //! \details This function calculates the attenuation coefficients
  //! from the attenuation image.
  //! \warning This is a temp function, probably a new class
  //! should be created, instead.
  Succeeded
    calculate_multiplicative_projdata(shared_ptr<ProjData> template_proj_data_ptr,
                           shared_ptr<VoxelsOnCartesianGrid<float> > atten_image_sptr,
                           shared_ptr<ProjData> out_proj_data_ptr,
                           std::string output_filename = "");


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

  //! remove cached attenuation integrals
  /*! should be used before recalculating scatter for a new attenuation image or
    when changing the sampling of the detector etc */
  virtual void remove_cache_for_integrals_over_attenuation();

  //! reset cached activity integrals
  /*! should be used before recalculating scatter for a new activity image or
    when changing the sampling of the detector etc */
  virtual void remove_cache_for_integrals_over_activity();

 private:

  //!
  //! \brief _scatter_simulation
  //! \details Class which will implement the scatter simulation.
  shared_ptr < ScatterSimulation > _scatter_simulation;

  Array<2,float> cached_activity_integral_scattpoint_det;
  Array<2,float> cached_attenuation_integral_scattpoint_det;


  //! set-up cache for attenuation integrals
  /*! \warning This will not remove existing cached data (if the sizes match). If you need this,
      call remove_cache_for_scattpoint_det_integrals_over_attenuation() first.
  */
  void initialise_cache_for_scattpoint_det_integrals_over_attenuation();
  //! set-up cache for activity integrals
  /*! \warning This will not remove existing cached data (if the sizes match). If you need this,
      call remove_cache_for_scattpoint_det_integrals_over_activity() first.
  */
  void initialise_cache_for_scattpoint_det_integrals_over_activity();

  //!
  //! \brief dets_to_voxs
  //! \param _num
  //! \param _axis
  //! \return
  //! \details A helper function which finds the optimal voxels size
  //! for a detector size
  inline float
  num_dets_to_vox_size(int _num, bool _axis);

  //!
  //! \brief voxs_to_dets
  //! \param _num
  //! \param _axis
  //! \return
  //! \details A helper function which finds the optimal det size
  //! for a voxel size.
  inline int
  vox_size_to_num_dets(float _num, bool _axis);


};


END_NAMESPACE_STIR

#include "stir/scatter/ScatterEstimationByBin.inl"

#endif
