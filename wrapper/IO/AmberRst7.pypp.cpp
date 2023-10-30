// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "AmberRst7.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/generalunitproperty.h"

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireBase/timeproperty.h"

#include "SireIO/amberformat.h"

#include "SireIO/amberrst7.h"

#include "SireIO/errors.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomvelocities.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidx.h"

#include "SireMol/trajectory.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "amberrst7.h"

#include "amberrst7.h"

SireIO::AmberRst7 __copy__(const SireIO::AmberRst7 &other){ return SireIO::AmberRst7(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AmberRst7_class(){

    { //::SireIO::AmberRst7
        typedef bp::class_< SireIO::AmberRst7, bp::bases< SireIO::MoleculeParser, SireBase::Property > > AmberRst7_exposer_t;
        AmberRst7_exposer_t AmberRst7_exposer = AmberRst7_exposer_t( "AmberRst7", "This class represents an Amber-format restartcoordinate file (ascii),\ncurrently supporting these files from Amber7 to Amber16.\n\nThe format of this file is described here;\n\nhttp:ambermd.orgformats.html\n\n(specifically the AMBER coordinaterestart file specification\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope AmberRst7_scope( AmberRst7_exposer );
        AmberRst7_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct by parsing the passed file") );
        AmberRst7_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct by parsing the data in the passed text lines") );
        AmberRst7_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct by extracting the necessary data from the passed System") );
        AmberRst7_exposer.def( bp::init< SireIO::AmberRst7 const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::AmberRst7::boxAngles
        
            typedef ::QVector< SireUnits::Dimension::PhysUnit< 0, 0, 0, 0, 0, 0, 1 > > ( ::SireIO::AmberRst7::*boxAngles_function_type)(  ) const;
            boxAngles_function_type boxAngles_function_value( &::SireIO::AmberRst7::boxAngles );
            
            AmberRst7_exposer.def( 
                "boxAngles"
                , boxAngles_function_value
                , bp::release_gil_policy()
                , "Return the parsed box angles" );
        
        }
        { //::SireIO::AmberRst7::boxDimensions
        
            typedef ::SireMaths::Vector ( ::SireIO::AmberRst7::*boxDimensions_function_type)(  ) const;
            boxDimensions_function_type boxDimensions_function_value( &::SireIO::AmberRst7::boxDimensions );
            
            AmberRst7_exposer.def( 
                "boxDimensions"
                , boxDimensions_function_value
                , bp::release_gil_policy()
                , "Return the parsed box dimensions" );
        
        }
        { //::SireIO::AmberRst7::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::AmberRst7::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::AmberRst7::construct );
            
            AmberRst7_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed filename" );
        
        }
        { //::SireIO::AmberRst7::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::AmberRst7::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::AmberRst7::construct );
            
            AmberRst7_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed set of lines" );
        
        }
        { //::SireIO::AmberRst7::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::AmberRst7::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::AmberRst7::construct );
            
            AmberRst7_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed SireSystem::System" );
        
        }
        { //::SireIO::AmberRst7::coordinates
        
            typedef ::QVector< SireMaths::Vector > ( ::SireIO::AmberRst7::*coordinates_function_type)(  ) const;
            coordinates_function_type coordinates_function_value( &::SireIO::AmberRst7::coordinates );
            
            AmberRst7_exposer.def( 
                "coordinates"
                , coordinates_function_value
                , bp::release_gil_policy()
                , "Return the parsed coordinate data" );
        
        }
        { //::SireIO::AmberRst7::formatDescription
        
            typedef ::QString ( ::SireIO::AmberRst7::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::AmberRst7::formatDescription );
            
            AmberRst7_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , bp::release_gil_policy()
                , "Return a description of the file format" );
        
        }
        { //::SireIO::AmberRst7::formatName
        
            typedef ::QString ( ::SireIO::AmberRst7::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::AmberRst7::formatName );
            
            AmberRst7_exposer.def( 
                "formatName"
                , formatName_function_value
                , bp::release_gil_policy()
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::AmberRst7::formatSuffix
        
            typedef ::QStringList ( ::SireIO::AmberRst7::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::AmberRst7::formatSuffix );
            
            AmberRst7_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , bp::release_gil_policy()
                , "Return the suffixes that RST7 files will typically have" );
        
        }
        { //::SireIO::AmberRst7::getFrame
        
            typedef ::SireMol::Frame ( ::SireIO::AmberRst7::*getFrame_function_type)( int ) const;
            getFrame_function_type getFrame_function_value( &::SireIO::AmberRst7::getFrame );
            
            AmberRst7_exposer.def( 
                "getFrame"
                , getFrame_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::AmberRst7::hasVelocities
        
            typedef bool ( ::SireIO::AmberRst7::*hasVelocities_function_type)(  ) const;
            hasVelocities_function_type hasVelocities_function_value( &::SireIO::AmberRst7::hasVelocities );
            
            AmberRst7_exposer.def( 
                "hasVelocities"
                , hasVelocities_function_value
                , bp::release_gil_policy()
                , "Return whether or not this restart file also provides velocities" );
        
        }
        { //::SireIO::AmberRst7::isFrame
        
            typedef bool ( ::SireIO::AmberRst7::*isFrame_function_type)(  ) const;
            isFrame_function_type isFrame_function_value( &::SireIO::AmberRst7::isFrame );
            
            AmberRst7_exposer.def( 
                "isFrame"
                , isFrame_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::AmberRst7::nAtoms
        
            typedef int ( ::SireIO::AmberRst7::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::AmberRst7::nAtoms );
            
            AmberRst7_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms whose coordinates are contained in this restart file" );
        
        }
        { //::SireIO::AmberRst7::nFrames
        
            typedef int ( ::SireIO::AmberRst7::*nFrames_function_type)(  ) const;
            nFrames_function_type nFrames_function_value( &::SireIO::AmberRst7::nFrames );
            
            AmberRst7_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AmberRst7_exposer.def( bp::self != bp::self );
        { //::SireIO::AmberRst7::operator=
        
            typedef ::SireIO::AmberRst7 & ( ::SireIO::AmberRst7::*assign_function_type)( ::SireIO::AmberRst7 const & ) ;
            assign_function_type assign_function_value( &::SireIO::AmberRst7::operator= );
            
            AmberRst7_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AmberRst7_exposer.def( bp::self == bp::self );
        { //::SireIO::AmberRst7::parse
        
            typedef ::SireIO::AmberRst7 ( *parse_function_type )( ::QString const & );
            parse_function_type parse_function_value( &::SireIO::AmberRst7::parse );
            
            AmberRst7_exposer.def( 
                "parse"
                , parse_function_value
                , ( bp::arg("filename") )
                , bp::release_gil_policy()
                , "Parse from the passed file" );
        
        }
        { //::SireIO::AmberRst7::time
        
            typedef double ( ::SireIO::AmberRst7::*time_function_type)(  ) const;
            time_function_type time_function_value( &::SireIO::AmberRst7::time );
            
            AmberRst7_exposer.def( 
                "time"
                , time_function_value
                , bp::release_gil_policy()
                , "Return the current time of the simulation from which this restart\nfile was written in picoseconds.\nThis is a negative number if the time has not been set" );
        
        }
        { //::SireIO::AmberRst7::title
        
            typedef ::QString ( ::SireIO::AmberRst7::*title_function_type)(  ) const;
            title_function_type title_function_value( &::SireIO::AmberRst7::title );
            
            AmberRst7_exposer.def( 
                "title"
                , title_function_value
                , bp::release_gil_policy()
                , "Return the title of the file" );
        
        }
        { //::SireIO::AmberRst7::toString
        
            typedef ::QString ( ::SireIO::AmberRst7::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::AmberRst7::toString );
            
            AmberRst7_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::AmberRst7::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::AmberRst7::typeName );
            
            AmberRst7_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::AmberRst7::velocities
        
            typedef ::QVector< SireMaths::Vector > ( ::SireIO::AmberRst7::*velocities_function_type)(  ) const;
            velocities_function_type velocities_function_value( &::SireIO::AmberRst7::velocities );
            
            AmberRst7_exposer.def( 
                "velocities"
                , velocities_function_value
                , bp::release_gil_policy()
                , "Return the parsed coordinate data" );
        
        }
        { //::SireIO::AmberRst7::what
        
            typedef char const * ( ::SireIO::AmberRst7::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::AmberRst7::what );
            
            AmberRst7_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AmberRst7_exposer.staticmethod( "parse" );
        AmberRst7_exposer.staticmethod( "typeName" );
        AmberRst7_exposer.def( "__copy__", &__copy__);
        AmberRst7_exposer.def( "__deepcopy__", &__copy__);
        AmberRst7_exposer.def( "clone", &__copy__);
        AmberRst7_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::AmberRst7 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AmberRst7_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::AmberRst7 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AmberRst7_exposer.def_pickle(sire_pickle_suite< ::SireIO::AmberRst7 >());
        AmberRst7_exposer.def( "__str__", &__str__< ::SireIO::AmberRst7 > );
        AmberRst7_exposer.def( "__repr__", &__str__< ::SireIO::AmberRst7 > );
    }

}
