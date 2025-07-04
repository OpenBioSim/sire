// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ConnectivityEditor.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/console.h"

#include "SireBase/errors.h"

#include "SireBase/parallel.h"

#include "SireError/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "angleid.h"

#include "atomidxmapping.h"

#include "atommatcher.h"

#include "atomselection.h"

#include "bondid.h"

#include "connectivity.h"

#include "dihedralid.h"

#include "improperid.h"

#include "moleculedata.h"

#include "moleculeinfo.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "tostring.h"

#include <QDataStream>

#include <QDebug>

#include <QElapsedTimer>

#include <boost/assert.hpp>

#include "connectivity.h"

SireMol::ConnectivityEditor __copy__(const SireMol::ConnectivityEditor &other){ return SireMol::ConnectivityEditor(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_ConnectivityEditor_class(){

    { //::SireMol::ConnectivityEditor
        typedef bp::class_< SireMol::ConnectivityEditor, bp::bases< SireMol::ConnectivityBase, SireMol::MolViewProperty, SireBase::Property > > ConnectivityEditor_exposer_t;
        ConnectivityEditor_exposer_t ConnectivityEditor_exposer = ConnectivityEditor_exposer_t( "ConnectivityEditor", "An editor that can be used to edit a Connectivity object\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope ConnectivityEditor_scope( ConnectivityEditor_exposer );
        ConnectivityEditor_exposer.def( bp::init< SireMol::Connectivity const & >(( bp::arg("connectivity") ), "Construct an editor to edit a copy of the passed\nConnectivity object") );
        ConnectivityEditor_exposer.def( bp::init< SireMol::ConnectivityEditor const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::ConnectivityEditor::commit
        
            typedef ::SireMol::Connectivity ( ::SireMol::ConnectivityEditor::*commit_function_type)(  ) const;
            commit_function_type commit_function_value( &::SireMol::ConnectivityEditor::commit );
            
            ConnectivityEditor_exposer.def( 
                "commit"
                , commit_function_value
                , bp::release_gil_policy()
                , "Return the editied connectivity" );
        
        }
        { //::SireMol::ConnectivityEditor::connect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*connect_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx ) ;
            connect_function_type connect_function_value( &::SireMol::ConnectivityEditor::connect );
            
            ConnectivityEditor_exposer.def( 
                "connect"
                , connect_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") )
                , bp::return_self< >()
                , "Record the connection between the atoms at indicies atom0\nand atom1\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::connect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*connect_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) ;
            connect_function_type connect_function_value( &::SireMol::ConnectivityEditor::connect );
            
            ConnectivityEditor_exposer.def( 
                "connect"
                , connect_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") )
                , bp::return_self< >()
                , "Record a connection between the atom identified by atom0 and\nthe atom identified by atom1\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::connect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*connect_function_type)( ::SireMol::BondID const & ) ;
            connect_function_type connect_function_value( &::SireMol::ConnectivityEditor::connect );
            
            ConnectivityEditor_exposer.def( 
                "connect"
                , connect_function_value
                , ( bp::arg("bond") )
                , bp::return_self< >()
                , "Create a connection for the passed bond" );
        
        }
        { //::SireMol::ConnectivityEditor::connect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*connect_function_type)( ::QList< SireMol::BondID > const & ) ;
            connect_function_type connect_function_value( &::SireMol::ConnectivityEditor::connect );
            
            ConnectivityEditor_exposer.def( 
                "connect"
                , connect_function_value
                , ( bp::arg("bonds") )
                , bp::return_self< >()
                , "Create a connection for the passed bonds" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnect_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx ) ;
            disconnect_function_type disconnect_function_value( &::SireMol::ConnectivityEditor::disconnect );
            
            ConnectivityEditor_exposer.def( 
                "disconnect"
                , disconnect_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") )
                , bp::return_self< >()
                , "Remove the connection between the atoms at indicies atom0\nand atom1 - this does nothing if there isnt already a connection\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnect_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) ;
            disconnect_function_type disconnect_function_value( &::SireMol::ConnectivityEditor::disconnect );
            
            ConnectivityEditor_exposer.def( 
                "disconnect"
                , disconnect_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") )
                , bp::return_self< >()
                , "Disconnect the atoms that are identified by atom0 and atom1 -\nthis does nothing if there isnt a connection between these atoms\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnect_function_type)( ::SireMol::BondID const & ) ;
            disconnect_function_type disconnect_function_value( &::SireMol::ConnectivityEditor::disconnect );
            
            ConnectivityEditor_exposer.def( 
                "disconnect"
                , disconnect_function_value
                , ( bp::arg("bond") )
                , bp::return_self< >()
                , "Disconnect the atoms in the passed bond - this does nothing if the\n  atoms arent connected" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnect_function_type)( ::QList< SireMol::BondID > const & ) ;
            disconnect_function_type disconnect_function_value( &::SireMol::ConnectivityEditor::disconnect );
            
            ConnectivityEditor_exposer.def( 
                "disconnect"
                , disconnect_function_value
                , ( bp::arg("bonds") )
                , bp::return_self< >()
                , "Disconnect the atoms in the passed bonds - this does nothing for any\n  of the atoms that arent connected" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnect
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnect_function_type)( ::QList< SireMol::AtomIdx > const &,bool ) ;
            disconnect_function_type disconnect_function_value( &::SireMol::ConnectivityEditor::disconnect );
            
            ConnectivityEditor_exposer.def( 
                "disconnect"
                , disconnect_function_value
                , ( bp::arg("atoms"), bp::arg("exclusive")=(bool)(true) )
                , bp::return_self< >()
                , "Disconnect any and all bonds involving the passed atoms. If exclusive is true,\n  then this only removes connection where both atoms are in atoms, otherwise\n  it removes connections which have one of more atoms in atoms\n" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnectAll
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnectAll_function_type)( ::SireMol::AtomIdx ) ;
            disconnectAll_function_type disconnectAll_function_value( &::SireMol::ConnectivityEditor::disconnectAll );
            
            ConnectivityEditor_exposer.def( 
                "disconnectAll"
                , disconnectAll_function_value
                , ( bp::arg("atomidx") )
                , bp::return_self< >()
                , "Remove all of the connections to the atom at index atomidx\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnectAll
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnectAll_function_type)( ::SireMol::ResIdx ) ;
            disconnectAll_function_type disconnectAll_function_value( &::SireMol::ConnectivityEditor::disconnectAll );
            
            ConnectivityEditor_exposer.def( 
                "disconnectAll"
                , disconnectAll_function_value
                , ( bp::arg("residx") )
                , bp::return_self< >()
                , "Remove all of the connections that involve any of the atoms\nin the residue at index residx\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnectAll
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnectAll_function_type)( ::SireMol::AtomID const & ) ;
            disconnectAll_function_type disconnectAll_function_value( &::SireMol::ConnectivityEditor::disconnectAll );
            
            ConnectivityEditor_exposer.def( 
                "disconnectAll"
                , disconnectAll_function_value
                , ( bp::arg("atomid") )
                , bp::return_self< >()
                , "Remove all of the connections to the atom identified by atomid\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnectAll
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnectAll_function_type)( ::SireMol::ResID const & ) ;
            disconnectAll_function_type disconnectAll_function_value( &::SireMol::ConnectivityEditor::disconnectAll );
            
            ConnectivityEditor_exposer.def( 
                "disconnectAll"
                , disconnectAll_function_value
                , ( bp::arg("resid") )
                , bp::return_self< >()
                , "Remove all of the connections that involve any of the atoms\nin the residue identified by resid\nThrow: SireMol::missing_residue\nThrow: SireMol::duplicate_residue\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ConnectivityEditor::disconnectAll
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*disconnectAll_function_type)(  ) ;
            disconnectAll_function_type disconnectAll_function_value( &::SireMol::ConnectivityEditor::disconnectAll );
            
            ConnectivityEditor_exposer.def( 
                "disconnectAll"
                , disconnectAll_function_value
                , bp::return_self< >()
                , "Remove all bonds from this molecule" );
        
        }
        ConnectivityEditor_exposer.def( bp::self != bp::self );
        { //::SireMol::ConnectivityEditor::operator=
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*assign_function_type)( ::SireMol::ConnectivityBase const & ) ;
            assign_function_type assign_function_value( &::SireMol::ConnectivityEditor::operator= );
            
            ConnectivityEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ConnectivityEditor_exposer.def( bp::self == bp::self );
        { //::SireMol::ConnectivityEditor::removeProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*removeProperty_function_type)( ::QString const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::ConnectivityEditor::removeProperty );
            
            ConnectivityEditor_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("key") )
                , bp::return_self< >()
                , "Remove the specified property from all bonds" );
        
        }
        { //::SireMol::ConnectivityEditor::removeProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*removeProperty_function_type)( ::SireMol::BondID const &,::QString const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::ConnectivityEditor::removeProperty );
            
            ConnectivityEditor_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("bond"), bp::arg("key") )
                , bp::return_self< >()
                , "Remove the specified property from the specified bond" );
        
        }
        { //::SireMol::ConnectivityEditor::removeProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*removeProperty_function_type)( ::SireMol::AngleID const &,::QString const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::ConnectivityEditor::removeProperty );
            
            ConnectivityEditor_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("ang"), bp::arg("key") )
                , bp::return_self< >()
                , "Remove the specified property from the specified angle" );
        
        }
        { //::SireMol::ConnectivityEditor::removeProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*removeProperty_function_type)( ::SireMol::DihedralID const &,::QString const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::ConnectivityEditor::removeProperty );
            
            ConnectivityEditor_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("dih"), bp::arg("key") )
                , bp::return_self< >()
                , "Remove the specified property from the specified dihedral" );
        
        }
        { //::SireMol::ConnectivityEditor::removeProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*removeProperty_function_type)( ::SireMol::ImproperID const &,::QString const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMol::ConnectivityEditor::removeProperty );
            
            ConnectivityEditor_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("imp"), bp::arg("key") )
                , bp::return_self< >()
                , "Remove the specified property from the specified improper" );
        
        }
        { //::SireMol::ConnectivityEditor::setProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*setProperty_function_type)( ::SireMol::BondID const &,::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMol::ConnectivityEditor::setProperty );
            
            ConnectivityEditor_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("bond"), bp::arg("key"), bp::arg("value") )
                , bp::return_self< >()
                , "Set the property for the specified bond, at the specified key, to value" );
        
        }
        { //::SireMol::ConnectivityEditor::setProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*setProperty_function_type)( ::SireMol::AngleID const &,::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMol::ConnectivityEditor::setProperty );
            
            ConnectivityEditor_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("ang"), bp::arg("key"), bp::arg("value") )
                , bp::return_self< >()
                , "Set the property for the specified angle, at the specified key, to value" );
        
        }
        { //::SireMol::ConnectivityEditor::setProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*setProperty_function_type)( ::SireMol::DihedralID const &,::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMol::ConnectivityEditor::setProperty );
            
            ConnectivityEditor_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("dih"), bp::arg("key"), bp::arg("value") )
                , bp::return_self< >()
                , "Set the property for the specified dihedral, at the specified key, to value" );
        
        }
        { //::SireMol::ConnectivityEditor::setProperty
        
            typedef ::SireMol::ConnectivityEditor & ( ::SireMol::ConnectivityEditor::*setProperty_function_type)( ::SireMol::ImproperID const &,::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMol::ConnectivityEditor::setProperty );
            
            ConnectivityEditor_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("imp"), bp::arg("key"), bp::arg("value") )
                , bp::return_self< >()
                , "Set the property for the specified improper, at the specified key, to value" );
        
        }
        { //::SireMol::ConnectivityEditor::takeProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::ConnectivityEditor::*takeProperty_function_type)( ::SireMol::BondID const &,::QString const & ) ;
            takeProperty_function_type takeProperty_function_value( &::SireMol::ConnectivityEditor::takeProperty );
            
            ConnectivityEditor_exposer.def( 
                "takeProperty"
                , takeProperty_function_value
                , ( bp::arg("bond"), bp::arg("key") )
                , bp::release_gil_policy()
                , "Take the specified property from the specified bond - this removes\nand returns the property if it exists. If it doesnt, then\na NullProperty is returned\n" );
        
        }
        { //::SireMol::ConnectivityEditor::takeProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::ConnectivityEditor::*takeProperty_function_type)( ::SireMol::AngleID const &,::QString const & ) ;
            takeProperty_function_type takeProperty_function_value( &::SireMol::ConnectivityEditor::takeProperty );
            
            ConnectivityEditor_exposer.def( 
                "takeProperty"
                , takeProperty_function_value
                , ( bp::arg("ang"), bp::arg("key") )
                , bp::release_gil_policy()
                , "Take the specified property from the specified angle - this removes\nand returns the property if it exists. If it doesnt, then\na NullProperty is returned\n" );
        
        }
        { //::SireMol::ConnectivityEditor::takeProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::ConnectivityEditor::*takeProperty_function_type)( ::SireMol::DihedralID const &,::QString const & ) ;
            takeProperty_function_type takeProperty_function_value( &::SireMol::ConnectivityEditor::takeProperty );
            
            ConnectivityEditor_exposer.def( 
                "takeProperty"
                , takeProperty_function_value
                , ( bp::arg("dih"), bp::arg("key") )
                , bp::release_gil_policy()
                , "Take the specified property from the specified dihedral - this removes\nand returns the property if it exists. If it doesnt, then\na NullProperty is returned\n" );
        
        }
        { //::SireMol::ConnectivityEditor::takeProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::ConnectivityEditor::*takeProperty_function_type)( ::SireMol::ImproperID const &,::QString const & ) ;
            takeProperty_function_type takeProperty_function_value( &::SireMol::ConnectivityEditor::takeProperty );
            
            ConnectivityEditor_exposer.def( 
                "takeProperty"
                , takeProperty_function_value
                , ( bp::arg("imp"), bp::arg("key") )
                , bp::release_gil_policy()
                , "Take the specified property from the specified improper - this removes\nand returns the property if it exists. If it doesnt, then\na NullProperty is returned\n" );
        
        }
        { //::SireMol::ConnectivityEditor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ConnectivityEditor::typeName );
            
            ConnectivityEditor_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ConnectivityEditor_exposer.staticmethod( "typeName" );
        ConnectivityEditor_exposer.def( "__copy__", &__copy__<SireMol::ConnectivityEditor>);
        ConnectivityEditor_exposer.def( "__deepcopy__", &__copy__<SireMol::ConnectivityEditor>);
        ConnectivityEditor_exposer.def( "clone", &__copy__<SireMol::ConnectivityEditor>);
        ConnectivityEditor_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ConnectivityEditor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ConnectivityEditor_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ConnectivityEditor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ConnectivityEditor_exposer.def_pickle(sire_pickle_suite< ::SireMol::ConnectivityEditor >());
        ConnectivityEditor_exposer.def( "__str__", &__str__< ::SireMol::ConnectivityEditor > );
        ConnectivityEditor_exposer.def( "__repr__", &__str__< ::SireMol::ConnectivityEditor > );
    }

}
