// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GroTop.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/booleanproperty.h"

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireIO/grotop.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "grotop.h"

#include <QFileInfo>

#include <QRegExp>

#include "grotop.h"

SireIO::GroTop __copy__(const SireIO::GroTop &other){ return SireIO::GroTop(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_GroTop_class(){

    { //::SireIO::GroTop
        typedef bp::class_< SireIO::GroTop, bp::bases< SireIO::MoleculeParser, SireBase::Property > > GroTop_exposer_t;
        GroTop_exposer_t GroTop_exposer = GroTop_exposer_t( "GroTop", "This class holds a parser for reading and writing Gromacs top topology files.\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope GroTop_scope( GroTop_exposer );
        GroTop_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the file called filename. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        GroTop_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the passed text lines. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        GroTop_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct this parser by extracting all necessary information from the\npassed SireSystem::System, looking for the properties that are specified\nin the passed property map") );
        GroTop_exposer.def( bp::init< SireIO::GroTop const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::GroTop::combiningRules
        
            typedef int ( ::SireIO::GroTop::*combiningRules_function_type)(  ) const;
            combiningRules_function_type combiningRules_function_value( &::SireIO::GroTop::combiningRules );
            
            GroTop_exposer.def( 
                "combiningRules"
                , combiningRules_function_value
                , "" );
        
        }
        { //::SireIO::GroTop::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::GroTop::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::GroTop::construct );
            
            GroTop_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\nfile using the passed properties" );
        
        }
        { //::SireIO::GroTop::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::GroTop::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::GroTop::construct );
            
            GroTop_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\ntext lines using the passed properties" );
        
        }
        { //::SireIO::GroTop::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::GroTop::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::GroTop::construct );
            
            GroTop_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , "Return the parser that has been constructed by extract all necessary\ndata from the passed SireSystem::System using the specified properties" );
        
        }
        { //::SireIO::GroTop::formatDescription
        
            typedef ::QString ( ::SireIO::GroTop::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::GroTop::formatDescription );
            
            GroTop_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , "Return a description of the file format" );
        
        }
        { //::SireIO::GroTop::formatName
        
            typedef ::QString ( ::SireIO::GroTop::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::GroTop::formatName );
            
            GroTop_exposer.def( 
                "formatName"
                , formatName_function_value
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::GroTop::formatSuffix
        
            typedef ::QStringList ( ::SireIO::GroTop::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::GroTop::formatSuffix );
            
            GroTop_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , "Return the suffixes that these files are normally associated with" );
        
        }
        { //::SireIO::GroTop::fudgeLJ
        
            typedef double ( ::SireIO::GroTop::*fudgeLJ_function_type)(  ) const;
            fudgeLJ_function_type fudgeLJ_function_value( &::SireIO::GroTop::fudgeLJ );
            
            GroTop_exposer.def( 
                "fudgeLJ"
                , fudgeLJ_function_value
                , "" );
        
        }
        { //::SireIO::GroTop::fudgeQQ
        
            typedef double ( ::SireIO::GroTop::*fudgeQQ_function_type)(  ) const;
            fudgeQQ_function_type fudgeQQ_function_value( &::SireIO::GroTop::fudgeQQ );
            
            GroTop_exposer.def( 
                "fudgeQQ"
                , fudgeQQ_function_value
                , "" );
        
        }
        { //::SireIO::GroTop::generateNonBondedPairs
        
            typedef bool ( ::SireIO::GroTop::*generateNonBondedPairs_function_type)(  ) const;
            generateNonBondedPairs_function_type generateNonBondedPairs_function_value( &::SireIO::GroTop::generateNonBondedPairs );
            
            GroTop_exposer.def( 
                "generateNonBondedPairs"
                , generateNonBondedPairs_function_value
                , "" );
        
        }
        { //::SireIO::GroTop::includePath
        
            typedef ::QStringList ( ::SireIO::GroTop::*includePath_function_type)( bool ) const;
            includePath_function_type includePath_function_value( &::SireIO::GroTop::includePath );
            
            GroTop_exposer.def( 
                "includePath"
                , includePath_function_value
                , ( bp::arg("absolute_paths")=(bool)(false) )
                , "Return the list of names of directories in which to search for\ninclude files. The directories are either absolute, or relative\nto the current directory. If absolute_paths is true then\nthe full absolute paths for directories that exist on this\nmachine will be returned" );
        
        }
        { //::SireIO::GroTop::includedFiles
        
            typedef ::QStringList ( ::SireIO::GroTop::*includedFiles_function_type)( bool ) const;
            includedFiles_function_type includedFiles_function_value( &::SireIO::GroTop::includedFiles );
            
            GroTop_exposer.def( 
                "includedFiles"
                , includedFiles_function_value
                , ( bp::arg("absolute_paths")=(bool)(false) )
                , "Return the list of names of files that were included when reading or\nwriting this file. The files are relative. If absolute_paths\nis true then the full absolute paths for the files will be\nused" );
        
        }
        { //::SireIO::GroTop::nonBondedFunctionType
        
            typedef int ( ::SireIO::GroTop::*nonBondedFunctionType_function_type)(  ) const;
            nonBondedFunctionType_function_type nonBondedFunctionType_function_value( &::SireIO::GroTop::nonBondedFunctionType );
            
            GroTop_exposer.def( 
                "nonBondedFunctionType"
                , nonBondedFunctionType_function_value
                , "" );
        
        }
        GroTop_exposer.def( bp::self != bp::self );
        { //::SireIO::GroTop::operator=
        
            typedef ::SireIO::GroTop & ( ::SireIO::GroTop::*assign_function_type)( ::SireIO::GroTop const & ) ;
            assign_function_type assign_function_value( &::SireIO::GroTop::operator= );
            
            GroTop_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GroTop_exposer.def( bp::self == bp::self );
        { //::SireIO::GroTop::toString
        
            typedef ::QString ( ::SireIO::GroTop::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::GroTop::toString );
            
            GroTop_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this parser" );
        
        }
        { //::SireIO::GroTop::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::GroTop::typeName );
            
            GroTop_exposer.def( 
                "typeName"
                , typeName_function_value
                , "Return the C++ name for this class" );
        
        }
        { //::SireIO::GroTop::what
        
            typedef char const * ( ::SireIO::GroTop::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::GroTop::what );
            
            GroTop_exposer.def( 
                "what"
                , what_function_value
                , "Return the C++ name for this class" );
        
        }
        GroTop_exposer.staticmethod( "typeName" );
        GroTop_exposer.def( "__copy__", &__copy__);
        GroTop_exposer.def( "__deepcopy__", &__copy__);
        GroTop_exposer.def( "clone", &__copy__);
        GroTop_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::GroTop >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GroTop_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::GroTop >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GroTop_exposer.def( "__str__", &__str__< ::SireIO::GroTop > );
        GroTop_exposer.def( "__repr__", &__str__< ::SireIO::GroTop > );
    }

}
