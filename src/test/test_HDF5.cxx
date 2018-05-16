/*
    Copyright (C) 2018, University of Hull
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
  \brief A simple program to test reading and orientation of GE Signa
  Normalisation factors.

  The GE toolbox was used to extract a set of nomarlisation factors, which
  are then read in Arrays.

  Then the BinNormalisationFromGEHDF5 is used to read the full normalisation
  data and these two are compared.

  * The geometric factors are stored in 3D arrays of size 16x357x1981. In this
    test only a 2x357x1981 is used.

  \author Nikos Efthimiou
*/
#include "stir/IndexRange.h"
#include "stir/IndexRange2D.h"
#include "stir/IndexRange3D.h"
#include "stir/Coordinate3D.h"
#include "stir/Scanner.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/recon_buildblock/BinNormalisationFromGEHDF5.h"
#include "stir/recon_buildblock/BinNormalisationFromProjData.h"
#include "stir/recon_buildblock/TrivialBinNormalisation.h"
#include "stir/ProjDataInMemory.h"
#include "stir/RunTests.h"
#include "stir/NumericType.h"
#include "stir/ByteOrder.h"
#include "stir/info.h"
#include "stir/warning.h"
#include "stir/IO/read_data.h"

using std::cerr;
using std::endl;

START_NAMESPACE_STIR

/*!
  \brief Class with tests for IndexRange, IndexRange3D.
*/
class HDF5_tests : public RunTests
{
public:
    void run_tests();
protected:
    //! Allocate arrays
    void allocate_space();
    //! Load the arrays with the normalisation factors for checking
    void load_check_norms();

    //! Array for the geometric factors.
    std::unique_ptr<Array<2,short> > check_geo;
    //! Template of the GE Signa.
    shared_ptr<Scanner>  ge_scanner;

    shared_ptr<ProjDataInfo> ge_proj_data_info;
    //! The number of tangestial positions
    int num_tang_poss;
    //! The number of axial / segment positions.
    int num_axial_poss;
    //! The number of view;
    //! In the test only 2 are used, in order to keep the sizes small.
    int num_view_poss;
};


void
HDF5_tests::run_tests()
{

    // Initialise the scanner
    ge_scanner.reset(new Scanner(Scanner::PETMR_Signa));

    // initialise a ProjDataInfo
    // TODO: span 2
    ge_proj_data_info.reset(
                dynamic_cast<ProjDataInfoCylindricalNoArcCorr *>(
                    ProjDataInfo::ProjDataInfoCTI(ge_scanner,
                                                  /*span=*/1, ge_scanner->get_num_rings()-1,
                                                  /*num_views,=*/ge_scanner->get_num_detectors_per_ring()/2,
                                                  /*num_tangential_poss=*/ge_scanner->get_max_num_non_arccorrected_bins(),
                                                  /*arc_corrected =*/false) ));

    // Exam Data
    shared_ptr<ExamInfo> exam_info_sptr(new ExamInfo());

    // Set sizes;
    num_axial_poss = 0;

    num_tang_poss = ge_proj_data_info->get_num_tangential_poss();

    for (int i = ge_proj_data_info->get_min_segment_num();
         i<= ge_proj_data_info->get_max_segment_num(); ++i)
        num_axial_poss += ge_proj_data_info->get_num_axial_poss(i);

    num_view_poss = 1; //! \todo: make this 2

    // Allocate space
    allocate_space();

    // Load defaults
    load_check_norms();

    // BinNormalization
    const std::string norm_fileName("/home/nikos/Desktop/norm3d");
    shared_ptr<BinNormalisation> bin_norm_sptr( new BinNormalisationFromGEHDF5(norm_fileName));


    ProjDataInMemory testNorm(exam_info_sptr, ge_proj_data_info, 1);
    testNorm.fill(1.f);

    // Apply or Undo
    bin_norm_sptr->apply(testNorm,0.0, 1.0);
    //bin_norm_sptr->undo(testNorm,0.0, 1.0);

    // Comparison between the check_geo which we have loaded from the disk.
    // and it has the correct orientation and values.
    // with the ProjData which have been corrected by the norm3d.

    //This will reveal the segments order. As you know which segements you have stored in
    // the check_geo.

    std::vector<int> segmentsOrder(ge_proj_data_info->get_num_segments());

    //   1.  [0] - [-1][1] - [-2][2]
    //    or
    //   2.  [0] - [1][-1]- [2][-2]

    // make the segmentOrder hold the values in the diagram above. Either 1. or 2.
    // will be correct.
    //    for (i=0; i < ge_proj_data_info->get_num_segments(); ++i)
    //    {
    //        segmentsOrder
    //    }

    for (int i_seg=0; i_seg < segmentsOrder.size(); ++i_seg)
    {
        for (int i_ax = ge_proj_data_info->get_min_axial_pos_num(i_seg);
            i_ax <= ge_proj_data_info->get_max_axial_pos_num(i_seg); ++i_ax)
        {
            for (int i_tang = ge_proj_data_info->get_min_tangential_pos_num();
                 i_tang <= ge_proj_data_info->get_max_tangential_pos_num(); ++i_tang)
            {
                // fill a  norm_geo Array. Initially will be a 2d array and then a 3D for the final version of
                // the test
            }
        }
    }

    // Finally compare the check_geo with the norm_geo. If they are the same ... cool.
    // Otherwise try the different segment order.

    // In case of a doubt try to print them in a file and compare visually.
 }

void
HDF5_tests::allocate_space()
{
    IndexRange2D ind(0, num_tang_poss-1, 0, 1980);
    check_geo.reset(new Array<2, short>(ind));
}

void
HDF5_tests::load_check_norms()
{
    std::string fileName("/home/nikos/Desktop/test_HDF5_nom_geo_dd.txt");
    std::ifstream inData(fileName, std::ifstream::in);

    if (read_data(inData, *check_geo) == Succeeded::no)
        error("Failed loading check geometric factors");

    //Print them for check. -- Commented out

    int nikos = 0;
}

END_NAMESPACE_STIR


USING_NAMESPACE_STIR
int main()
{
    HDF5_tests tests;
    tests.run_tests();
    return tests.main_return_value();

}
