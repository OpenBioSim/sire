// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "IntraGroupFF.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/booleanproperty.h"

#include "SireBase/errors.h"

#include "SireBase/lengthproperty.h"

#include "SireBase/refcountdata.h"

#include "SireError/errors.h"

#include "SireMol/atomselection.h"

#include "SireMol/molecule.h"

#include "SireMol/molecules.h"

#include "SireMol/molresid.h"

#include "SireMol/partialmolecule.h"

#include "SireMol/residue.h"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "cljcalculator.h"

#include "cljshiftfunction.h"

#include "intragroupff.h"

#include <QDebug>

#include <QElapsedTimer>

#include "intragroupff.h"

SireMM::IntraGroupFF __copy__(const SireMM::IntraGroupFF &other){ return SireMM::IntraGroupFF(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_IntraGroupFF_class(){

    { //::SireMM::IntraGroupFF
        typedef bp::class_< SireMM::IntraGroupFF, bp::bases< SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > IntraGroupFF_exposer_t;
        IntraGroupFF_exposer_t IntraGroupFF_exposer = IntraGroupFF_exposer_t( "IntraGroupFF", "This forcefield is used to calculate the intramolecular\ncoulomb and LJ energy of the contained molecules. Note\nthat this is the coulomb and LJ energy of the non-bonded\natoms only, i.e. it does not contain the scaled\n1-4 coulomb and LJ energies. These should be calculated\nseparately, e.g. via additional terms added to InternalFF\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope IntraGroupFF_scope( IntraGroupFF_exposer );
        IntraGroupFF_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "Construct, specifying the name of the forcefield") );
        IntraGroupFF_exposer.def( bp::init< SireMM::IntraGroupFF const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::IntraGroupFF::accept
        
            typedef void ( ::SireMM::IntraGroupFF::*accept_function_type)(  ) ;
            accept_function_type accept_function_value( &::SireMM::IntraGroupFF::accept );
            
            IntraGroupFF_exposer.def( 
                "accept"
                , accept_function_value
                , bp::release_gil_policy()
                , "Tell the forcefield that the last move was accepted. This tells the\nforcefield to make permanent any temporary changes that were used a workspace\nto avoid memory allocation during a move" );
        
        }
        { //::SireMM::IntraGroupFF::cljFunction
        
            typedef ::SireMM::CLJIntraFunction const & ( ::SireMM::IntraGroupFF::*cljFunction_function_type)(  ) const;
            cljFunction_function_type cljFunction_function_value( &::SireMM::IntraGroupFF::cljFunction );
            
            IntraGroupFF_exposer.def( 
                "cljFunction"
                , cljFunction_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the function used to calculate the energy" );
        
        }
        { //::SireMM::IntraGroupFF::cljFunction
        
            typedef ::SireMM::CLJIntraFunction const & ( ::SireMM::IntraGroupFF::*cljFunction_function_type)( ::QString ) const;
            cljFunction_function_type cljFunction_function_value( &::SireMM::IntraGroupFF::cljFunction );
            
            IntraGroupFF_exposer.def( 
                "cljFunction"
                , cljFunction_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the CLJFunction associated with the passed key" );
        
        }
        { //::SireMM::IntraGroupFF::cljFunctionKeys
        
            typedef ::QStringList ( ::SireMM::IntraGroupFF::*cljFunctionKeys_function_type)(  ) const;
            cljFunctionKeys_function_type cljFunctionKeys_function_value( &::SireMM::IntraGroupFF::cljFunctionKeys );
            
            IntraGroupFF_exposer.def( 
                "cljFunctionKeys"
                , cljFunctionKeys_function_value
                , bp::release_gil_policy()
                , "Return the keys of all CLJFunctions added to this forcefield" );
        
        }
        { //::SireMM::IntraGroupFF::cljFunctions
        
            typedef ::QHash< QString, SireBase::PropPtr< SireMM::CLJFunction > > ( ::SireMM::IntraGroupFF::*cljFunctions_function_type)(  ) const;
            cljFunctions_function_type cljFunctions_function_value( &::SireMM::IntraGroupFF::cljFunctions );
            
            IntraGroupFF_exposer.def( 
                "cljFunctions"
                , cljFunctions_function_value
                , bp::release_gil_policy()
                , "Return the hash of all CLJFunctions in this forcefield, indexed by their key" );
        
        }
        { //::SireMM::IntraGroupFF::components
        
            typedef ::SireMM::MultiCLJComponent const & ( ::SireMM::IntraGroupFF::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireMM::IntraGroupFF::components );
            
            IntraGroupFF_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the energy components of this forcefield" );
        
        }
        { //::SireMM::IntraGroupFF::containsProperty
        
            typedef bool ( ::SireMM::IntraGroupFF::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::IntraGroupFF::containsProperty );
            
            IntraGroupFF_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Return whether or not this forcefield contains the property property" );
        
        }
        { //::SireMM::IntraGroupFF::disableParallelCalculation
        
            typedef void ( ::SireMM::IntraGroupFF::*disableParallelCalculation_function_type)(  ) ;
            disableParallelCalculation_function_type disableParallelCalculation_function_value( &::SireMM::IntraGroupFF::disableParallelCalculation );
            
            IntraGroupFF_exposer.def( 
                "disableParallelCalculation"
                , disableParallelCalculation_function_value
                , bp::release_gil_policy()
                , "Turn off use of a multicore parallel calculation of the energy.\nThis may be quicker if you have few atoms in the forcefield,\nor if you are only planning on allocating one core per forcefield" );
        
        }
        { //::SireMM::IntraGroupFF::disableReproducibleCalculation
        
            typedef void ( ::SireMM::IntraGroupFF::*disableReproducibleCalculation_function_type)(  ) ;
            disableReproducibleCalculation_function_type disableReproducibleCalculation_function_value( &::SireMM::IntraGroupFF::disableReproducibleCalculation );
            
            IntraGroupFF_exposer.def( 
                "disableReproducibleCalculation"
                , disableReproducibleCalculation_function_value
                , bp::release_gil_policy()
                , "Turn off an energy summing algorithm that guarantees the same energy\nregardless of whether a single core or multicore calculation is being\nperformed (i.e. rounding errors in both cases will not be identical)" );
        
        }
        { //::SireMM::IntraGroupFF::enableParallelCalculation
        
            typedef void ( ::SireMM::IntraGroupFF::*enableParallelCalculation_function_type)(  ) ;
            enableParallelCalculation_function_type enableParallelCalculation_function_value( &::SireMM::IntraGroupFF::enableParallelCalculation );
            
            IntraGroupFF_exposer.def( 
                "enableParallelCalculation"
                , enableParallelCalculation_function_value
                , bp::release_gil_policy()
                , "Turn on use of a multicore parallel calculation of the energy.\nThis is on by default, and spreads the energy calculations over\navailable cores" );
        
        }
        { //::SireMM::IntraGroupFF::enableReproducibleCalculation
        
            typedef void ( ::SireMM::IntraGroupFF::*enableReproducibleCalculation_function_type)(  ) ;
            enableReproducibleCalculation_function_type enableReproducibleCalculation_function_value( &::SireMM::IntraGroupFF::enableReproducibleCalculation );
            
            IntraGroupFF_exposer.def( 
                "enableReproducibleCalculation"
                , enableReproducibleCalculation_function_value
                , bp::release_gil_policy()
                , "Turn on an energy summing algorithm that guarantees the same energy\nregardless of whether a single core or multicore calculation is being\nperformed (i.e. rounding errors in both cases will be identical)" );
        
        }
        { //::SireMM::IntraGroupFF::mustNowRecalculateFromScratch
        
            typedef void ( ::SireMM::IntraGroupFF::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireMM::IntraGroupFF::mustNowRecalculateFromScratch );
            
            IntraGroupFF_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , bp::release_gil_policy()
                , "Signal that this forcefield must now be recalculated from scratch" );
        
        }
        { //::SireMM::IntraGroupFF::nCLJFunctions
        
            typedef int ( ::SireMM::IntraGroupFF::*nCLJFunctions_function_type)(  ) const;
            nCLJFunctions_function_type nCLJFunctions_function_value( &::SireMM::IntraGroupFF::nCLJFunctions );
            
            IntraGroupFF_exposer.def( 
                "nCLJFunctions"
                , nCLJFunctions_function_value
                , bp::release_gil_policy()
                , "Return the number of CLJ functions in this forcefield. There should always\nbe at least one" );
        
        }
        { //::SireMM::IntraGroupFF::needsAccepting
        
            typedef bool ( ::SireMM::IntraGroupFF::*needsAccepting_function_type)(  ) const;
            needsAccepting_function_type needsAccepting_function_value( &::SireMM::IntraGroupFF::needsAccepting );
            
            IntraGroupFF_exposer.def( 
                "needsAccepting"
                , needsAccepting_function_value
                , bp::release_gil_policy()
                , "Return whether or not this forcefield is using a temporary workspace that\nneeds to be accepted" );
        
        }
        IntraGroupFF_exposer.def( bp::self != bp::self );
        { //::SireMM::IntraGroupFF::operator=
        
            typedef ::SireMM::IntraGroupFF & ( ::SireMM::IntraGroupFF::*assign_function_type)( ::SireMM::IntraGroupFF const & ) ;
            assign_function_type assign_function_value( &::SireMM::IntraGroupFF::operator= );
            
            IntraGroupFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IntraGroupFF_exposer.def( bp::self == bp::self );
        { //::SireMM::IntraGroupFF::properties
        
            typedef ::SireBase::Properties const & ( ::SireMM::IntraGroupFF::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::IntraGroupFF::properties );
            
            IntraGroupFF_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return all of the properties of this function" );
        
        }
        { //::SireMM::IntraGroupFF::property
        
            typedef ::SireBase::Property const & ( ::SireMM::IntraGroupFF::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::IntraGroupFF::property );
            
            IntraGroupFF_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the value of the forcefield property with name name" );
        
        }
        { //::SireMM::IntraGroupFF::removeAllCLJFunctions
        
            typedef void ( ::SireMM::IntraGroupFF::*removeAllCLJFunctions_function_type)(  ) ;
            removeAllCLJFunctions_function_type removeAllCLJFunctions_function_value( &::SireMM::IntraGroupFF::removeAllCLJFunctions );
            
            IntraGroupFF_exposer.def( 
                "removeAllCLJFunctions"
                , removeAllCLJFunctions_function_value
                , bp::release_gil_policy()
                , "Function to remove all of the CLJFunctions (except for the default function)" );
        
        }
        { //::SireMM::IntraGroupFF::removeCLJFunctionAt
        
            typedef void ( ::SireMM::IntraGroupFF::*removeCLJFunctionAt_function_type)( ::QString ) ;
            removeCLJFunctionAt_function_type removeCLJFunctionAt_function_value( &::SireMM::IntraGroupFF::removeCLJFunctionAt );
            
            IntraGroupFF_exposer.def( 
                "removeCLJFunctionAt"
                , removeCLJFunctionAt_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Remove the CLJ function with key key - note that you cannot remove\nthe default CLJ function" );
        
        }
        { //::SireMM::IntraGroupFF::setCLJFunction
        
            typedef void ( ::SireMM::IntraGroupFF::*setCLJFunction_function_type)( ::SireMM::CLJIntraFunction const & ) ;
            setCLJFunction_function_type setCLJFunction_function_value( &::SireMM::IntraGroupFF::setCLJFunction );
            
            IntraGroupFF_exposer.def( 
                "setCLJFunction"
                , setCLJFunction_function_value
                , ( bp::arg("cljfunc") )
                , bp::release_gil_policy()
                , "Function used to set the CLJIntraFunction used to calculate\nthe intramolecular energy" );
        
        }
        { //::SireMM::IntraGroupFF::setCLJFunction
        
            typedef void ( ::SireMM::IntraGroupFF::*setCLJFunction_function_type)( ::QString,::SireMM::CLJIntraFunction const & ) ;
            setCLJFunction_function_type setCLJFunction_function_value( &::SireMM::IntraGroupFF::setCLJFunction );
            
            IntraGroupFF_exposer.def( 
                "setCLJFunction"
                , setCLJFunction_function_value
                , ( bp::arg("key"), bp::arg("cljfunc") )
                , bp::release_gil_policy()
                , "Set the CLJFunction with key key equal to cljfunc" );
        
        }
        { //::SireMM::IntraGroupFF::setProperty
        
            typedef bool ( ::SireMM::IntraGroupFF::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMM::IntraGroupFF::setProperty );
            
            IntraGroupFF_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") )
                , bp::release_gil_policy()
                , "Set the forcefield property called name to the value property. Note that\nthis only affects the default CLJFunction. Additional functions must\nbe configured before adding them to the forcefield" );
        
        }
        { //::SireMM::IntraGroupFF::setUseParallelCalculation
        
            typedef void ( ::SireMM::IntraGroupFF::*setUseParallelCalculation_function_type)( bool ) ;
            setUseParallelCalculation_function_type setUseParallelCalculation_function_value( &::SireMM::IntraGroupFF::setUseParallelCalculation );
            
            IntraGroupFF_exposer.def( 
                "setUseParallelCalculation"
                , setUseParallelCalculation_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Set whether or not to use a multicore parallel algorithm\nto calculate the energy" );
        
        }
        { //::SireMM::IntraGroupFF::setUseReproducibleCalculation
        
            typedef void ( ::SireMM::IntraGroupFF::*setUseReproducibleCalculation_function_type)( bool ) ;
            setUseReproducibleCalculation_function_type setUseReproducibleCalculation_function_value( &::SireMM::IntraGroupFF::setUseReproducibleCalculation );
            
            IntraGroupFF_exposer.def( 
                "setUseReproducibleCalculation"
                , setUseReproducibleCalculation_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Switch on or off use of an energy summing algorithm that guarantees the\nsame energy regardless of whether a single core or multicore calculation\nis being performed" );
        
        }
        { //::SireMM::IntraGroupFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::IntraGroupFF::typeName );
            
            IntraGroupFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::IntraGroupFF::usesParallelCalculation
        
            typedef bool ( ::SireMM::IntraGroupFF::*usesParallelCalculation_function_type)(  ) const;
            usesParallelCalculation_function_type usesParallelCalculation_function_value( &::SireMM::IntraGroupFF::usesParallelCalculation );
            
            IntraGroupFF_exposer.def( 
                "usesParallelCalculation"
                , usesParallelCalculation_function_value
                , bp::release_gil_policy()
                , "Return whether or not a parallel algorithm is used to calculate energies" );
        
        }
        { //::SireMM::IntraGroupFF::usesReproducibleCalculation
        
            typedef bool ( ::SireMM::IntraGroupFF::*usesReproducibleCalculation_function_type)(  ) const;
            usesReproducibleCalculation_function_type usesReproducibleCalculation_function_value( &::SireMM::IntraGroupFF::usesReproducibleCalculation );
            
            IntraGroupFF_exposer.def( 
                "usesReproducibleCalculation"
                , usesReproducibleCalculation_function_value
                , bp::release_gil_policy()
                , "Return whether or not a reproducible energy summing algorithm is being\nused to accumulate the energies" );
        
        }
        { //::SireMM::IntraGroupFF::what
        
            typedef char const * ( ::SireMM::IntraGroupFF::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::IntraGroupFF::what );
            
            IntraGroupFF_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IntraGroupFF_exposer.staticmethod( "typeName" );
        IntraGroupFF_exposer.def( "__copy__", &__copy__<SireMM::IntraGroupFF>);
        IntraGroupFF_exposer.def( "__deepcopy__", &__copy__<SireMM::IntraGroupFF>);
        IntraGroupFF_exposer.def( "clone", &__copy__<SireMM::IntraGroupFF>);
        IntraGroupFF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::IntraGroupFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupFF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::IntraGroupFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupFF_exposer.def_pickle(sire_pickle_suite< ::SireMM::IntraGroupFF >());
        IntraGroupFF_exposer.def( "__str__", &__str__< ::SireMM::IntraGroupFF > );
        IntraGroupFF_exposer.def( "__repr__", &__str__< ::SireMM::IntraGroupFF > );
        IntraGroupFF_exposer.def( "__len__", &__len_count< ::SireMM::IntraGroupFF > );
    }

}
