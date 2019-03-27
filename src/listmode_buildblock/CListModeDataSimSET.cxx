/*
    Copyright (C) 2019 University of Hull
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
  \ingroup listmode
  \brief Implementation of class stir::CListModeDataSimSET

  \author Nikos Efthimiou
*/

#include "stir/listmode/CListModeDataSimSET.h"
#include "stir/Scanner.h"
#include "stir/Succeeded.h"
#include "stir/FilePath.h"
#include "stir/info.h"
#include "stir/warning.h"
#include "stir/error.h"
#include <boost/format.hpp>

extern "C" {
#include <LbEnvironment.h>
#include <LbInterface.h>
}

#define	PHGRDHST_NumFlags	3
/* LOCAL CONSTANTS */
#define PHGRDHST_IsUsePHGHistory()		LbFgIsSet(PhgOptions, LBFlag0)		/* Will we use the PHG history file */
#define PHGRDHST_IsUseColHistory()		LbFgIsSet(PhgOptions, LBFlag1)		/* Will we use the Collimator history file */
#define PHGRDHST_IsUseDetHistory()		LbFgIsSet(PhgOptions, LBFlag2)		/* Will we use the Detector history file */

START_NAMESPACE_STIR

//! N.E: Large parts adapted from phgbin and functions called by it;
CListModeDataSimSET::
CListModeDataSimSET(const std::string& _phg_filename)
    : phg_filename(_phg_filename)
{
    set_defaults();

    double SubObjCurTimeBinDuration = 0.0;

    // The dirty trick we did in SimSETListmodeInputFileFormat::is_SimSET_signature()
    char** argv = new char*[3];
    argv[0] = nullptr;
    argv[1] = nullptr;
    argv[2] = nullptr;

    char* pseudo_binary = new char[7];
    memset(pseudo_binary, 0, 7);
    strcpy(pseudo_binary, "phgbin\0");
    argv[0] = pseudo_binary;

    char* flag = new char[4];
    memset(flag, 0, 4);
    strcpy(flag, "-p");
    flag[3] = '\0';
    argv[1] = flag;

    char *argv_c = new char[phg_filename.size() + 1];
    memset(argv_c, 0, phg_filename.size()+1);
    phg_filename.copy(argv_c, phg_filename.size());
    argv_c[phg_filename.size()] = '\0';
    argv[2] = argv_c;

    char *knownOptions = new char[4];
    strcpy(knownOptions, "pcd");
    knownOptions[3] = '\0';

    char optArgs[PHGRDHST_NumFlags][LBEnMxArgLen];
    LbUsFourByte optArgFlags = (LBFlag0);
    LbUsFourByte phgrdhstArgIndex;

    /* Perform initialization tasks */

    /* Get our options */
    if (!LbEnGetOptions(3, argv, &knownOptions,
                        &PhgOptions, optArgs, optArgFlags, &phgrdhstArgIndex))
    {
        error("CListModeDataSimSET: Unable to get options.");
    }

    /* Make sure the didn't specify more than one history file */
    if (PHGRDHST_IsUsePHGHistory() && (PHGRDHST_IsUseColHistory() || PHGRDHST_IsUseDetHistory()))
    {
        error("CListModeDataSimSET: You can only specify one type of history file.");
    }

    /* Make sure they specified a history file */
    // N.E.: This should never be the case
    if (!PHGRDHST_IsUsePHGHistory() && !PHGRDHST_IsUseColHistory() && !PHGRDHST_IsUseDetHistory())
    {
        error("CListModeDataSimSET: You must specify the use of a PHG history file (-p)\n"
              " or a collimator history file (-c)\n"
              " or a detector history file (-d)\n");
    }

    strcpy(PhgRunTimeParams.PhgParamFilePath,argv[2]);

    /* Get our run-time parameters */
    if (!PhgGetRunTimeParams())
        error("CListModeDataSimSET: Error geting our run-time parameters.");

    if (PhgRunTimeParams.PhgIsPET != 1)
        error("CListModeDataSimSET: Currently we are able to process only PET files.");

    strcpy(phgrdhstHistName, PhgRunTimeParams.PhgPhoHFileHistoryFilePath);
    strcpy(phgrdhstHistParamsName, PhgRunTimeParams.PhgPhoHParamsFilePath);

    phgrdhstHistName_str.push_back(*phgrdhstHistName);
    phgrdhstHistParamsName_str.push_back(*phgrdhstHistParamsName);

    /* Initialize the math library - NE: skipped*/
    /* Initialize the emission list manager - NE: skipped*/
    /* Initialize the sub-object manager - NE: partial*/

    /* Set the length of the current time bin */
    SubObjCurTimeBinDuration = static_cast<double>(PhgRunTimeParams.Phg_LengthOfScan);

    /* Initialize the productivity table manager - NE: skipped*/
    /* Initialize the Cylinder Positions */
    {
        /* Set object cylinder */
        if (!CylPosInitObjectCylinder())
        {
            error("CListModeDataSimSET: Error initialize the Cylinder Positions.");
        }

        /* Set the criticial zone */
        if (!CylPosInitCriticalZone(PhgRunTimeParams.Phg_AcceptanceAngle))
        {
            error("CListModeDataSimSET: Error setting the criticial zone.");
        }

        /* Set the limit cylinder */
        CylPosInitLimitCylinder();
    }

    /* Setup up the productivity information - NE: Reluctantly skipped */
    /* Initialize the photon tracking module - NE: skipped */
    // In the following part several nonPET and simulation parts, are skipped.
    /* We support only PET files so this should leed to an error always. */
    ColCurParams = 0;
    if (PHG_IsCollimateOnTheFly())
    {
        error("CListModeDataSimSET: This history files has collimator information. \n"
              "Currently we support only PET simulations.");
    }

    /* Initialize the detection module if necessary */
    if (PHG_IsDetectOnTheFly())
    {
        if (!DetInitialize(PHGRDHST_IsUsePHGHistory()))
            error("CListModeDataSimSET: Unable to initialize the detection module.");
    }


    // Get Binning parameter file.
    //    char *binParams_cstr = phgrdhstHdrParams->H.PhgRunTimeParams.PhgBinParamsFilePath[PhgNumBinParams];

    // Try to guess scanner
    // CylPosTargetCylinder holds the goemetric information on the scanner.
    // //   double			radius;		/* Radius of cylinder */
    // //	double			zMin;		/* Minimum z coordinate of cylinder */
    // //	double			zMax;		/* Maximum z coordinate of cylinder */
    // //	double			centerX;	/* Center point on x axis */
    // //	double			centerY;	/* Center point on y axis */

    // ExamInfo initialisation
    this->exam_info_sptr.reset(new ExamInfo);

    // Only PET scanners supported
    this->exam_info_sptr->imaging_modality = ImagingModality::PT;
    this->exam_info_sptr->originating_system = this->originating_system;
    this->exam_info_sptr->set_low_energy_thres(static_cast<float>(PhgRunTimeParams.PhgMinimumEnergy));
//    this->exam_info_sptr->set_high_energy_thres(static_cast<float>(PhgRunTimeParams.PhgMinimumEnergy));


    // initialise ProjData.


    int nikos = 0;


    int crap = 0;

    // if the ProjData have been initialised properly create a
    // Input Stream from SimSET.
    history_file_sptr.reset(new InputStreamFromSimSET());
    history_file_sptr->set_up(phgrdhstHistParamsName_str);

    // Clean up.
    delete [] argv;
    delete [] pseudo_binary;
    delete [] argv_c;
    delete [] flag;
    delete [] knownOptions;
}

std::string
CListModeDataSimSET::
get_name() const
{
    return phg_filename;
}

CListModeDataSimSET::~CListModeDataSimSET()
{
    delete [] phgrdhstHistName;
    delete [] phgrdhstHistParamsName;
}

shared_ptr <CListRecord>
CListModeDataSimSET::
get_empty_record_sptr() const
{
    shared_ptr<CListRecord> sptr(new CListRecordSimSET(this->get_proj_data_info_sptr()->get_scanner_sptr()));
    return sptr;
}

Succeeded
CListModeDataSimSET::
open_lm_file()
{
    //    info(boost::format("CListModeDataSimSET: used ROOT file %s") %
    //         this->root_file_sptr->get_SimSET_filename());
    return Succeeded::yes;
}



Succeeded
CListModeDataSimSET::
get_next_record(CListRecord& record_of_general_type) const
{
    //    CListModeDataSimSET& record = dynamic_cast<CListModeDataSimSET&>(record_of_general_type);
    //    return root_file_sptr->get_next_record(record);
}

Succeeded
CListModeDataSimSET::
reset()
{
    return history_file_sptr->reset();
}

unsigned long CListModeDataSimSET::get_total_number_of_events() const
{
    //    return root_file_sptr->get_total_number_of_events();
}

CListModeData::SavedPosition
CListModeDataSimSET::
save_get_position()
{
    return static_cast<SavedPosition>(history_file_sptr->save_get_position());
}

Succeeded
CListModeDataSimSET::
set_get_position(const CListModeDataSimSET::SavedPosition& pos)
{
    return history_file_sptr->set_get_position(pos);
}

void
CListModeDataSimSET::
set_defaults()
{
    phgrdhstHistParamsName = new char[1024];
    phgrdhstHistName = new char[1024];

    /* Clear the file name parameters */
    phgrdhstHistParamsName[0] = '\0';
    phgrdhstHistName[0] = '\0';
}

Succeeded
CListModeDataSimSET::
check_scanner_match_geometry(std::string& ret, const shared_ptr<Scanner>& scanner_sptr)
{
    //    std::ostringstream stream("CListModeDataSimSET: The Scanner does not match the GATE geometry. Check: ");
    //    bool ok = true;

    //    if (scanner_sptr->get_num_rings() != root_file_sptr->get_num_rings())
    //    {
    //        stream << "the number of rings, ";
    //        ok = false;
    //    }

    //    if (scanner_sptr->get_num_detectors_per_ring() != root_file_sptr->get_num_dets_per_ring())
    //    {
    //        stream << "the number of detector per ring, ";
    //        ok = false;
    //    }

    //    if (scanner_sptr->get_num_axial_blocks_per_bucket() != root_file_sptr->get_num_axial_blocks_per_bucket_v())
    //    {
    //        stream << "the number of axial blocks per bucket, ";
    //        ok = false;
    //    }

    //    if(scanner_sptr->get_num_transaxial_blocks_per_bucket() != root_file_sptr->get_num_transaxial_blocks_per_bucket_v())
    //    {
    //        stream << "the number of transaxial blocks per bucket, ";
    //        ok = false;
    //    }

    //    if(scanner_sptr->get_num_axial_crystals_per_block() != root_file_sptr->get_num_axial_crystals_per_block_v())
    //    {
    //        stream << "the number of axial crystals per block, ";
    //        ok = false;
    //    }

    //    if(scanner_sptr->get_num_transaxial_crystals_per_block() != root_file_sptr->get_num_transaxial_crystals_per_block_v())
    //    {
    //        stream << "the number of transaxial crystals per block, ";
    //        ok = false;
    //    }

    //    if(scanner_sptr->get_num_axial_crystals_per_singles_unit() != root_file_sptr->get_num_axial_crystals_per_singles_unit())
    //    {
    //        stream << "the number of axial crystals per singles unit, ";
    //        ok = false;
    //    }

    //    if(scanner_sptr->get_num_transaxial_crystals_per_singles_unit() != root_file_sptr->get_num_trans_crystals_per_singles_unit())
    //    {
    //        stream << "the number of transaxial crystals per singles unit, ";
    //        ok = false;
    //    }

    //    if (!ok)
    //    {
    //        ret = stream.str();
    //        return Succeeded::no;
    //    }

    return Succeeded::yes;
}

Succeeded
CListModeDataSimSET::
check_scanner_definition(std::string& ret)
{
    //    if ( num_rings == -1 ||
    //         num_detectors_per_ring == -1 ||
    //         max_num_non_arccorrected_bins == -1 ||
    //         inner_ring_diameter == -1.f ||
    //         average_depth_of_interaction == -1.f ||
    //         ring_spacing == -.1f ||
    //         bin_size == -1.f )
    //    {
    //       std::ostringstream stream("CListModeDataSimSET: The User_defined_scanner has not been fully described.\nPlease include in the hroot:\n");

    //       if (num_rings == -1)
    //           stream << "Number of rings := \n";

    //       if (num_detectors_per_ring == -1)
    //           stream << "Number of detectors per ring := \n";

    //       if (max_num_non_arccorrected_bins == -1)
    //           stream << "Maximum number of non-arc-corrected bins := \n";

    //       if (inner_ring_diameter == -1)
    //           stream << "Inner ring diameter (cm) := \n";

    //       if (average_depth_of_interaction == -1)
    //           stream << "Average depth of interaction (cm) := \n";

    //       if (ring_spacing == -1)
    //           stream << "Distance between rings (cm) := \n";

    //       if (bin_size == -1)
    //           stream << "Default bin size (cm) := \n";

    //       ret = stream.str();

    //       return Succeeded::no;
    //    }

    return Succeeded::yes;
}


END_NAMESPACE_STIR
