/*
 
*/
/*!
  \file
  \ingroup listmode
  \brief Declaration of class stir::CListModeDataSimSET

  \author Nikos Efthimiou
  \author Harry Tsoumpas
  \author Kris Thielemans
*/

#ifndef __stir_listmode_CListModeDataSimSET_H__
#define __stir_listmode_CListModeDataSimSET_H__

#include "stir/listmode/CListModeData.h"
#include "stir/listmode/CListRecordSimSET.h"
#include "stir/IO/InputStreamFromSimSET.h"
#include "stir/KeyParser.h"

START_NAMESPACE_STIR


class CListModeDataSimSET : public CListModeData
{
public:
    //! construct from the filename of the Interfile header
    CListModeDataSimSET(const std::string& hroot_filename_prefix);

    //! returns the header filename
    virtual std::string
    get_name() const;
    //! Set private members default values;
    void set_defaults();

    virtual
    shared_ptr <CListRecord> get_empty_record_sptr() const;

    virtual
    Succeeded get_next_record(CListRecord& record) const;

    virtual
    Succeeded reset();

    virtual
    SavedPosition save_get_position();

    virtual
    Succeeded set_get_position(const SavedPosition&);

    virtual
    bool has_delayeds() const { return true; }

    virtual inline
    unsigned long int
    get_total_number_of_events() const ;

private:
    //! Check if the hroot contains a full scanner description
    Succeeded check_scanner_definition(std::string& ret);
    //! Check if the scanner_sptr matches the geometry in root_file_sptr
    Succeeded check_scanner_match_geometry(std::string& ret, const shared_ptr<Scanner>& scanner_sptr);

    //! The header file
    std::string hroot_filename;

    //! Pointer to the listmode data
    shared_ptr<InputStreamFromSimSET > root_file_sptr;

//! \name Variables that can be set in the hroot file to define a scanner's geometry.
//! They are compared to the Scanner  (if set)  and the InputStreamFromROOTFile
//! geometry, as given by the repeaters. Can be used to check for inconsistencies.
//@{
    //! The name of the originating scanner
    std::string originating_system;
    //! Number of rings, set in the hroot file (optional)
    int num_rings;
    //! Number of detectors per ring, set in the hroot file (optional)
    int num_detectors_per_ring;
    //! Number of non arc corrected bins, set in the hroot file (optional)
    int max_num_non_arccorrected_bins;
    //! Inner ring diameter, set in the hroot file (optional)
    float inner_ring_diameter;
    //! Average depth of interaction, set in the hroot file (optional)
    float average_depth_of_interaction;
    //! Ring spacing, set in the hroot file (optional)
    float ring_spacing;
    //! Bin size, set in the hroot file (optional)
    float bin_size;
//@}

    KeyParser parser;

    Succeeded open_lm_file();
};


END_NAMESPACE_STIR

#endif
