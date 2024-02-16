// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "EMLECallback.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireVol/triclinicbox.h"

#include "emle.h"

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireVol/triclinicbox.h"

#include "emle.h"

SireOpenMM::EMLECallback __copy__(const SireOpenMM::EMLECallback &other){ return SireOpenMM::EMLECallback(other); }

const char* pvt_get_name(const SireOpenMM::EMLECallback&){ return "SireOpenMM::EMLECallback";}

#include "Helpers/release_gil_policy.hpp"

void register_EMLECallback_class(){

    { //::SireOpenMM::EMLECallback
        typedef bp::class_< SireOpenMM::EMLECallback > EMLECallback_exposer_t;
        EMLECallback_exposer_t EMLECallback_exposer = EMLECallback_exposer_t( "EMLECallback", "A callback wrapper class to allow use of electrostatic embedding of\nmachine learning potentials via emle-engine.", bp::init< >("Default constructor.") );
        bp::scope EMLECallback_scope( EMLECallback_exposer );
        EMLECallback_exposer.def( bp::init< bp::api::object, bp::optional< QString > >(( bp::arg("arg0"), bp::arg("callback")="_sire_callback" ), "Constructor\nPar:am py_object\nA Python object that contains the callback function.\n\nPar:am callback\nThe name of a callback method that take the following arguments:\n- numbers_qm: A list of atomic numbers for the atoms in the ML region.\n- charges_mm: A list of the MM charges in mod electron charge.\n- xyz_qm: A vector of positions for the atoms in the ML region in Angstrom.\n- xyz_mm: A vector of positions for the atoms in the MM region in Angstrom.\n") );
        { //::SireOpenMM::EMLECallback::call
        
            typedef ::boost::tuples::tuple< double, QVector< QVector< double > >, QVector< QVector< double > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireOpenMM::EMLECallback::*call_function_type)( ::QVector< int >,::QVector< double >,::QVector< QVector< double > >,::QVector< QVector< double > > ) const;
            call_function_type call_function_value( &::SireOpenMM::EMLECallback::call );
            
            EMLECallback_exposer.def( 
                "call"
                , call_function_value
                , ( bp::arg("numbers_qm"), bp::arg("charges_mm"), bp::arg("xyz_qm"), bp::arg("xyz_mm") )
                , bp::release_gil_policy()
                , "Call the callback function.\nPar:am numbers_qm\nA vector of atomic numbers for the atoms in the ML region.\n\nPar:am charges_mm\nA vector of the charges on the MM atoms in mod electron charge.\n\nPar:am xyz_qm\nA vector of positions for the atoms in the ML region in Angstrom.\n\nPar:am xyz_mm\nA vector of positions for the atoms in the MM region in Angstrom.\n\nReturn:s\nA tuple containing:\n- The energy in kJmol.\n- A vector of forces for the QM atoms in kJmolnm.\n- A vector of forces for the MM atoms in kJmolnm.\n" );
        
        }
        { //::SireOpenMM::EMLECallback::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireOpenMM::EMLECallback::typeName );
            
            EMLECallback_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class." );
        
        }
        { //::SireOpenMM::EMLECallback::what
        
            typedef char const * ( ::SireOpenMM::EMLECallback::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireOpenMM::EMLECallback::what );
            
            EMLECallback_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class." );
        
        }
        EMLECallback_exposer.staticmethod( "typeName" );
        EMLECallback_exposer.def( "__copy__", &__copy__);
        EMLECallback_exposer.def( "__deepcopy__", &__copy__);
        EMLECallback_exposer.def( "clone", &__copy__);
        EMLECallback_exposer.def( "__str__", &pvt_get_name);
        EMLECallback_exposer.def( "__repr__", &pvt_get_name);
    }

}