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
#ifndef __stir_IO_SimSETListmodeInputFileFormat_h__
#define __stir_IO_SimSETListmodeInputFileFormat_h__


#include "stir/IO/InputFileFormat.h"
#include "stir/listmode/CListModeDataSimSET.h"
#include "stir/error.h"
#include <iostream>

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

//!
//! \brief The SimSETListmodeInputFileFormat class
//! \details Class to check the file signature of SiSET history files.
//! These are listmode files.
//!
//! The phg parameter file should be used here.
//! However, STIR passes only 1024 bytes. The solution is that the
//! first line of your phg file should always be the fileName of that
//! file. That way we can use SimSET functions to check for
//! consistency.
//!
//!
//! Three things are checked Header size, header kind,
//! and version of the header.
//!
//! \warning Althought the first two work as supposed, in the third
//! we need to skip several byte to reach the proper value.
//! In the original function I cannot find a similar requirement.
//!
//! \warning In addition, that float number with the header version
//! has to be read as big endian. Further investigation on this should be
//! carried out.
//!
//! \author Nikos Efthimiou
//!
class SimSETListmodeInputFileFormat :
        public InputFileFormat<CListModeData >
{
public:
    virtual const std::string
    get_name() const
    {
        return "SimSET_Histrory_file";
    }

    virtual unique_ptr<data_type>
    read_from_file(std::istream& input) const
    {
        error("read_from_file for SimSET listmode data with istream not implemented %s:%s. Sorry",
              __FILE__, __LINE__);
        return unique_ptr<data_type>();
    }

    virtual unique_ptr<data_type>
    read_from_file(const std::string& filename) const
    {
        return unique_ptr<data_type>(new CListModeDataSimSET(filename));
    }

protected:

    virtual
    bool actual_can_read(const FileSignature& signature,
                         std::istream& input) const
    {
        return this->is_SimSET_signature(input);
    }

    bool is_SimSET_signature(std::istream& signature) const
    {
        std::string argv_str;
        std::getline(signature, argv_str);
        char fileName[1024];

        // Remove everything prior to := from the line.
        argv_str = argv_str.substr(argv_str.find("=") + 1);
        strcpy(fileName, argv_str.c_str());

        // Prepend the -p for SimSET to know the type of file.
        argv_str.insert(0, "-p ");

        char *argv = new char[argv_str.length() + 1];
        strcpy(argv, argv_str.c_str());

        char *knownOptions[] = {"pcd"};
        char optArgs[PHGRDHST_NumFlags][LBEnMxArgLen];
        LbUsFourByte optArgFlags = (LBFlag0);
        LbUsFourByte phgrdhstArgIndex;

        /* Get our options */
        if (!LbEnGetOptions(2, &argv, knownOptions,
                            &PhgOptions, optArgs, optArgFlags, &phgrdhstArgIndex)) {

            delete [] argv;
            return false;
        }

        /* Make sure the didn't specify more than one history file */
        if (PHGRDHST_IsUsePHGHistory() && (PHGRDHST_IsUseColHistory() || PHGRDHST_IsUseDetHistory()))
        {
            warning("SimSETListmodeInputFileFormat: You can only specify one type of history file.");
            return false;
        }

        delete [] argv;
        return true;
    }
};

END_NAMESPACE_STIR

#endif
