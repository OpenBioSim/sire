// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "PDBx.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireIO/pdbx.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/core.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/trajectory.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "pdbx.h"

#include "pdbx.h"

SireIO::PDBx __copy__(const SireIO::PDBx &other){ return SireIO::PDBx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_PDBx_class(){

    { //::SireIO::PDBx
        typedef bp::class_< SireIO::PDBx, bp::bases< SireIO::MoleculeParser, SireBase::Property > > PDBx_exposer_t;
        PDBx_exposer_t PDBx_exposer = PDBx_exposer_t( "PDBx", "This class holds a parser for reading and writing PDBxmmcif files", bp::init< >("Constructor") );
        bp::scope PDBx_scope( PDBx_exposer );
        PDBx_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the file called filename. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        PDBx_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the passed text lines. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        PDBx_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct this parser by extracting all necessary information from the\npassed SireSystem::System, looking for the properties that are specified\nin the passed property map") );
        PDBx_exposer.def( bp::init< SireIO::PDBx const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::PDBx::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::PDBx::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::PDBx::construct );
            
            PDBx_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return the parser that has been constructed by reading in the passed\nfile using the passed properties" );
        
        }
        { //::SireIO::PDBx::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::PDBx::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::PDBx::construct );
            
            PDBx_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return the parser that has been constructed by reading in the passed\ntext lines using the passed properties" );
        
        }
        { //::SireIO::PDBx::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::PDBx::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::PDBx::construct );
            
            PDBx_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return the parser that has been constructed by extract all necessary\ndata from the passed SireSystem::System using the specified properties" );
        
        }
        { //::SireIO::PDBx::formatDescription
        
            typedef ::QString ( ::SireIO::PDBx::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::PDBx::formatDescription );
            
            PDBx_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , bp::release_gil_policy()
                , "Return a description of the file format" );
        
        }
        { //::SireIO::PDBx::formatName
        
            typedef ::QString ( ::SireIO::PDBx::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::PDBx::formatName );
            
            PDBx_exposer.def( 
                "formatName"
                , formatName_function_value
                , bp::release_gil_policy()
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::PDBx::formatSuffix
        
            typedef ::QStringList ( ::SireIO::PDBx::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::PDBx::formatSuffix );
            
            PDBx_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , bp::release_gil_policy()
                , "Return the suffixes that these files are normally associated with" );
        
        }
        { //::SireIO::PDBx::getFrame
        
            typedef ::SireMol::Frame ( ::SireIO::PDBx::*getFrame_function_type)( int ) const;
            getFrame_function_type getFrame_function_value( &::SireIO::PDBx::getFrame );
            
            PDBx_exposer.def( 
                "getFrame"
                , getFrame_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::PDBx::isFrame
        
            typedef bool ( ::SireIO::PDBx::*isFrame_function_type)(  ) const;
            isFrame_function_type isFrame_function_value( &::SireIO::PDBx::isFrame );
            
            PDBx_exposer.def( 
                "isFrame"
                , isFrame_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::PDBx::isTopology
        
            typedef bool ( ::SireIO::PDBx::*isTopology_function_type)(  ) const;
            isTopology_function_type isTopology_function_value( &::SireIO::PDBx::isTopology );
            
            PDBx_exposer.def( 
                "isTopology"
                , isTopology_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::PDBx::nAtoms
        
            typedef int ( ::SireIO::PDBx::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::PDBx::nAtoms );
            
            PDBx_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the total number of atoms in all molecules." );
        
        }
        { //::SireIO::PDBx::nFrames
        
            typedef int ( ::SireIO::PDBx::*nFrames_function_type)(  ) const;
            nFrames_function_type nFrames_function_value( &::SireIO::PDBx::nFrames );
            
            PDBx_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        PDBx_exposer.def( bp::self != bp::self );
        { //::SireIO::PDBx::operator=
        
            typedef ::SireIO::PDBx & ( ::SireIO::PDBx::*assign_function_type)( ::SireIO::PDBx const & ) ;
            assign_function_type assign_function_value( &::SireIO::PDBx::operator= );
            
            PDBx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PDBx_exposer.def( bp::self == bp::self );
        { //::SireIO::PDBx::toLines
        
            typedef ::QVector< QString > ( ::SireIO::PDBx::*toLines_function_type)(  ) const;
            toLines_function_type toLines_function_value( &::SireIO::PDBx::toLines );
            
            PDBx_exposer.def( 
                "toLines"
                , toLines_function_value
                , bp::release_gil_policy()
                , "Convert the the parsed data to a collection of PDBx record lines." );
        
        }
        { //::SireIO::PDBx::toString
        
            typedef ::QString ( ::SireIO::PDBx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::PDBx::toString );
            
            PDBx_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this parser" );
        
        }
        { //::SireIO::PDBx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::PDBx::typeName );
            
            PDBx_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class" );
        
        }
        { //::SireIO::PDBx::what
        
            typedef char const * ( ::SireIO::PDBx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::PDBx::what );
            
            PDBx_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class" );
        
        }
        PDBx_exposer.staticmethod( "typeName" );
        PDBx_exposer.def( "__copy__", &__copy__);
        PDBx_exposer.def( "__deepcopy__", &__copy__);
        PDBx_exposer.def( "clone", &__copy__);
        PDBx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::PDBx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PDBx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::PDBx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PDBx_exposer.def_pickle(sire_pickle_suite< ::SireIO::PDBx >());
        PDBx_exposer.def( "__str__", &__str__< ::SireIO::PDBx > );
        PDBx_exposer.def( "__repr__", &__str__< ::SireIO::PDBx > );
    }

}
