// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AmberParams.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireCAS/expression.h"

#include "SireCAS/sum.h"

#include "SireCAS/symbol.h"

#include "SireCAS/trigfuncs.h"

#include "SireCAS/values.h"

#include "SireError/errors.h"

#include "SireMM/cljnbpairs.h"

#include "SireMM/fouratomfunctions.h"

#include "SireMM/threeatomfunctions.h"

#include "SireMM/twoatomfunctions.h"

#include "SireMol/angleid.h"

#include "SireMol/atomidx.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/dihedralid.h"

#include "SireMol/improperid.h"

#include "SireMol/molecule.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "amberparams.h"

#include "amberparams.h"

SireMM::AmberParams __copy__(const SireMM::AmberParams &other){ return SireMM::AmberParams(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AmberParams_class(){

    { //::SireMM::AmberParams
        typedef bp::class_< SireMM::AmberParams, bp::bases< SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > AmberParams_exposer_t;
        AmberParams_exposer_t AmberParams_exposer = AmberParams_exposer_t( "AmberParams", "This class stores AMBER bonded force field parameters for\na collection of bonds, angles, dihedrals, impropers\nand 1-4 scaling factors.\n\nAuthor: Julien Michel  Christopher Woods\n", bp::init< >("Null Constructor") );
        bp::scope AmberParams_scope( AmberParams_exposer );
        AmberParams_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ), "Constructor for the passed molecule") );
        AmberParams_exposer.def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") ), "Constructor for the passed molecule") );
        AmberParams_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "Constructor for the passed molecule") );
        AmberParams_exposer.def( bp::init< SireMM::AmberParams const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::AmberParams::add
        
            typedef void ( ::SireMM::AmberParams::*add_function_type)( ::SireMol::AtomID const &,::SireUnits::Dimension::Charge,::SireUnits::Dimension::MolarMass,::SireMol::Element const &,::SireMM::LJParameter const &,::QString const &,::SireUnits::Dimension::Length,double,::QString const & ) ;
            add_function_type add_function_value( &::SireMM::AmberParams::add );
            
            AmberParams_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("atom"), bp::arg("charge"), bp::arg("mass"), bp::arg("element"), bp::arg("ljparam"), bp::arg("amber_type"), bp::arg("born_radius"), bp::arg("screening_parameter"), bp::arg("treechain") )
                , "Set the atom parameters for the specified atom to the provided values" );
        
        }
        { //::SireMM::AmberParams::add
        
            typedef void ( ::SireMM::AmberParams::*add_function_type)( ::SireMol::BondID const &,double,double,bool ) ;
            add_function_type add_function_value( &::SireMM::AmberParams::add );
            
            AmberParams_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("bond"), bp::arg("k"), bp::arg("r0"), bp::arg("includes_hydrogen") )
                , "" );
        
        }
        { //::SireMM::AmberParams::add
        
            typedef void ( ::SireMM::AmberParams::*add_function_type)( ::SireMol::AngleID const &,double,double,bool ) ;
            add_function_type add_function_value( &::SireMM::AmberParams::add );
            
            AmberParams_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("angle"), bp::arg("k"), bp::arg("theta0"), bp::arg("includes_hydrogen") )
                , "" );
        
        }
        { //::SireMM::AmberParams::add
        
            typedef void ( ::SireMM::AmberParams::*add_function_type)( ::SireMol::DihedralID const &,double,double,double,bool ) ;
            add_function_type add_function_value( &::SireMM::AmberParams::add );
            
            AmberParams_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("dihedral"), bp::arg("k"), bp::arg("periodicity"), bp::arg("phase"), bp::arg("includes_hydrogen") )
                , "" );
        
        }
        { //::SireMM::AmberParams::add
        
            typedef void ( ::SireMM::AmberParams::*add_function_type)( ::SireMol::ImproperID const &,double,double,double,bool ) ;
            add_function_type add_function_value( &::SireMM::AmberParams::add );
            
            AmberParams_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("improper"), bp::arg("v"), bp::arg("periodicity"), bp::arg("phase"), bp::arg("includes_hydrogen") )
                , "" );
        
        }
        { //::SireMM::AmberParams::addNB14
        
            typedef void ( ::SireMM::AmberParams::*addNB14_function_type)( ::SireMol::BondID const &,double,double ) ;
            addNB14_function_type addNB14_function_value( &::SireMM::AmberParams::addNB14 );
            
            AmberParams_exposer.def( 
                "addNB14"
                , addNB14_function_value
                , ( bp::arg("pair"), bp::arg("cscl"), bp::arg("ljscl") )
                , "" );
        
        }
        { //::SireMM::AmberParams::amberTypes
        
            typedef ::SireMol::AtomStringProperty ( ::SireMM::AmberParams::*amberTypes_function_type)(  ) const;
            amberTypes_function_type amberTypes_function_value( &::SireMM::AmberParams::amberTypes );
            
            AmberParams_exposer.def( 
                "amberTypes"
                , amberTypes_function_value
                , "Return all of the amber atom types" );
        
        }
        { //::SireMM::AmberParams::angleFunctions
        
            typedef ::SireMM::ThreeAtomFunctions ( ::SireMM::AmberParams::*angleFunctions_function_type)(  ) const;
            angleFunctions_function_type angleFunctions_function_value( &::SireMM::AmberParams::angleFunctions );
            
            AmberParams_exposer.def( 
                "angleFunctions"
                , angleFunctions_function_value
                , "Return all of the angle parameters converted to a set of ThreeAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::angleFunctions
        
            typedef ::SireMM::ThreeAtomFunctions ( ::SireMM::AmberParams::*angleFunctions_function_type)( ::SireCAS::Symbol const & ) const;
            angleFunctions_function_type angleFunctions_function_value( &::SireMM::AmberParams::angleFunctions );
            
            AmberParams_exposer.def( 
                "angleFunctions"
                , angleFunctions_function_value
                , ( bp::arg("THETA") )
                , "Return all of the angle parameters converted to a set of ThreeAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::angles
        
            typedef ::QHash< SireMol::AngleID, QPair< SireMM::AmberAngle, bool > > ( ::SireMM::AmberParams::*angles_function_type)(  ) const;
            angles_function_type angles_function_value( &::SireMM::AmberParams::angles );
            
            AmberParams_exposer.def( 
                "angles"
                , angles_function_value
                , "" );
        
        }
        { //::SireMM::AmberParams::bondFunctions
        
            typedef ::SireMM::TwoAtomFunctions ( ::SireMM::AmberParams::*bondFunctions_function_type)(  ) const;
            bondFunctions_function_type bondFunctions_function_value( &::SireMM::AmberParams::bondFunctions );
            
            AmberParams_exposer.def( 
                "bondFunctions"
                , bondFunctions_function_value
                , "Return all of the bond parameters converted to a set of TwoAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::bondFunctions
        
            typedef ::SireMM::TwoAtomFunctions ( ::SireMM::AmberParams::*bondFunctions_function_type)( ::SireCAS::Symbol const & ) const;
            bondFunctions_function_type bondFunctions_function_value( &::SireMM::AmberParams::bondFunctions );
            
            AmberParams_exposer.def( 
                "bondFunctions"
                , bondFunctions_function_value
                , ( bp::arg("R") )
                , "Return all of the bond parameters converted to a set of TwoAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::bonds
        
            typedef ::QHash< SireMol::BondID, QPair< SireMM::AmberBond, bool > > ( ::SireMM::AmberParams::*bonds_function_type)(  ) const;
            bonds_function_type bonds_function_value( &::SireMM::AmberParams::bonds );
            
            AmberParams_exposer.def( 
                "bonds"
                , bonds_function_value
                , "" );
        
        }
        { //::SireMM::AmberParams::charges
        
            typedef ::SireMol::AtomCharges ( ::SireMM::AmberParams::*charges_function_type)(  ) const;
            charges_function_type charges_function_value( &::SireMM::AmberParams::charges );
            
            AmberParams_exposer.def( 
                "charges"
                , charges_function_value
                , "Return the charges on the atoms" );
        
        }
        { //::SireMM::AmberParams::cljScaleFactors
        
            typedef ::SireMM::CLJNBPairs ( ::SireMM::AmberParams::*cljScaleFactors_function_type)(  ) const;
            cljScaleFactors_function_type cljScaleFactors_function_value( &::SireMM::AmberParams::cljScaleFactors );
            
            AmberParams_exposer.def( 
                "cljScaleFactors"
                , cljScaleFactors_function_value
                , "Return the CLJ nonbonded 1-4 scale factors for the molecule" );
        
        }
        { //::SireMM::AmberParams::connectivity
        
            typedef ::SireMol::Connectivity ( ::SireMM::AmberParams::*connectivity_function_type)(  ) const;
            connectivity_function_type connectivity_function_value( &::SireMM::AmberParams::connectivity );
            
            AmberParams_exposer.def( 
                "connectivity"
                , connectivity_function_value
                , "Return the connectivity of the molecule implied by the\nthe bonds" );
        
        }
        { //::SireMM::AmberParams::convert
        
            typedef ::SireMol::BondID ( ::SireMM::AmberParams::*convert_function_type)( ::SireMol::BondID const & ) const;
            convert_function_type convert_function_value( &::SireMM::AmberParams::convert );
            
            AmberParams_exposer.def( 
                "convert"
                , convert_function_value
                , ( bp::arg("bond") )
                , "Convert the passed BondID into AtomIdx IDs, sorted in index order" );
        
        }
        { //::SireMM::AmberParams::convert
        
            typedef ::SireMol::AngleID ( ::SireMM::AmberParams::*convert_function_type)( ::SireMol::AngleID const & ) const;
            convert_function_type convert_function_value( &::SireMM::AmberParams::convert );
            
            AmberParams_exposer.def( 
                "convert"
                , convert_function_value
                , ( bp::arg("angle") )
                , "Convert the passed AngleID into AtomIdx IDs, sorted in index order" );
        
        }
        { //::SireMM::AmberParams::convert
        
            typedef ::SireMol::DihedralID ( ::SireMM::AmberParams::*convert_function_type)( ::SireMol::DihedralID const & ) const;
            convert_function_type convert_function_value( &::SireMM::AmberParams::convert );
            
            AmberParams_exposer.def( 
                "convert"
                , convert_function_value
                , ( bp::arg("dihedral") )
                , "Convert the passed DihedralID into AtomIdx IDs, sorted in index order" );
        
        }
        { //::SireMM::AmberParams::convert
        
            typedef ::SireMol::ImproperID ( ::SireMM::AmberParams::*convert_function_type)( ::SireMol::ImproperID const & ) const;
            convert_function_type convert_function_value( &::SireMM::AmberParams::convert );
            
            AmberParams_exposer.def( 
                "convert"
                , convert_function_value
                , ( bp::arg("improper") )
                , "Convert the passed ImproperID into AtomIdx IDs, sorted in index order" );
        
        }
        { //::SireMM::AmberParams::dihedralFunctions
        
            typedef ::SireMM::FourAtomFunctions ( ::SireMM::AmberParams::*dihedralFunctions_function_type)(  ) const;
            dihedralFunctions_function_type dihedralFunctions_function_value( &::SireMM::AmberParams::dihedralFunctions );
            
            AmberParams_exposer.def( 
                "dihedralFunctions"
                , dihedralFunctions_function_value
                , "Return all of the dihedral parameters converted to a set of FourAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::dihedralFunctions
        
            typedef ::SireMM::FourAtomFunctions ( ::SireMM::AmberParams::*dihedralFunctions_function_type)( ::SireCAS::Symbol const & ) const;
            dihedralFunctions_function_type dihedralFunctions_function_value( &::SireMM::AmberParams::dihedralFunctions );
            
            AmberParams_exposer.def( 
                "dihedralFunctions"
                , dihedralFunctions_function_value
                , ( bp::arg("PHI") )
                , "Return all of the dihedral parameters converted to a set of FourAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::dihedrals
        
            typedef ::QHash< SireMol::DihedralID, QPair< SireMM::AmberDihedral, bool > > ( ::SireMM::AmberParams::*dihedrals_function_type)(  ) const;
            dihedrals_function_type dihedrals_function_value( &::SireMM::AmberParams::dihedrals );
            
            AmberParams_exposer.def( 
                "dihedrals"
                , dihedrals_function_value
                , "" );
        
        }
        { //::SireMM::AmberParams::elements
        
            typedef ::SireMol::AtomElements ( ::SireMM::AmberParams::*elements_function_type)(  ) const;
            elements_function_type elements_function_value( &::SireMM::AmberParams::elements );
            
            AmberParams_exposer.def( 
                "elements"
                , elements_function_value
                , "Return the atom elements" );
        
        }
        { //::SireMM::AmberParams::excludedAtoms
        
            typedef ::SireMM::CLJNBPairs ( ::SireMM::AmberParams::*excludedAtoms_function_type)(  ) const;
            excludedAtoms_function_type excludedAtoms_function_value( &::SireMM::AmberParams::excludedAtoms );
            
            AmberParams_exposer.def( 
                "excludedAtoms"
                , excludedAtoms_function_value
                , "Return the excluded atoms of the molecule. The returned\nobject has a matrix of all atom pairs, where the value\nis 0 for atom0-atom1 pairs that are to be excluded,\nand 1 for atom0-atom1 pairs that are to be included\nin the nonbonded calculation" );
        
        }
        { //::SireMM::AmberParams::gbRadii
        
            typedef ::SireMol::AtomRadii ( ::SireMM::AmberParams::*gbRadii_function_type)(  ) const;
            gbRadii_function_type gbRadii_function_value( &::SireMM::AmberParams::gbRadii );
            
            AmberParams_exposer.def( 
                "gbRadii"
                , gbRadii_function_value
                , "Return all of the Born radii of the atoms" );
        
        }
        { //::SireMM::AmberParams::gbScreening
        
            typedef ::SireMol::AtomFloatProperty ( ::SireMM::AmberParams::*gbScreening_function_type)(  ) const;
            gbScreening_function_type gbScreening_function_value( &::SireMM::AmberParams::gbScreening );
            
            AmberParams_exposer.def( 
                "gbScreening"
                , gbScreening_function_value
                , "Return all of the Born screening parameters for the atoms" );
        
        }
        { //::SireMM::AmberParams::getNB14
        
            typedef ::SireMM::AmberNB14 ( ::SireMM::AmberParams::*getNB14_function_type)( ::SireMol::BondID const & ) const;
            getNB14_function_type getNB14_function_value( &::SireMM::AmberParams::getNB14 );
            
            AmberParams_exposer.def( 
                "getNB14"
                , getNB14_function_value
                , ( bp::arg("pair") )
                , "" );
        
        }
        { //::SireMM::AmberParams::getParameter
        
            typedef ::SireMM::AmberBond ( ::SireMM::AmberParams::*getParameter_function_type)( ::SireMol::BondID const & ) const;
            getParameter_function_type getParameter_function_value( &::SireMM::AmberParams::getParameter );
            
            AmberParams_exposer.def( 
                "getParameter"
                , getParameter_function_value
                , ( bp::arg("bond") )
                , "" );
        
        }
        { //::SireMM::AmberParams::getParameter
        
            typedef ::SireMM::AmberAngle ( ::SireMM::AmberParams::*getParameter_function_type)( ::SireMol::AngleID const & ) const;
            getParameter_function_type getParameter_function_value( &::SireMM::AmberParams::getParameter );
            
            AmberParams_exposer.def( 
                "getParameter"
                , getParameter_function_value
                , ( bp::arg("angle") )
                , "" );
        
        }
        { //::SireMM::AmberParams::getParameter
        
            typedef ::SireMM::AmberDihedral ( ::SireMM::AmberParams::*getParameter_function_type)( ::SireMol::DihedralID const & ) const;
            getParameter_function_type getParameter_function_value( &::SireMM::AmberParams::getParameter );
            
            AmberParams_exposer.def( 
                "getParameter"
                , getParameter_function_value
                , ( bp::arg("dihedral") )
                , "" );
        
        }
        { //::SireMM::AmberParams::getParameter
        
            typedef ::SireMM::AmberDihedral ( ::SireMM::AmberParams::*getParameter_function_type)( ::SireMol::ImproperID const & ) const;
            getParameter_function_type getParameter_function_value( &::SireMM::AmberParams::getParameter );
            
            AmberParams_exposer.def( 
                "getParameter"
                , getParameter_function_value
                , ( bp::arg("improper") )
                , "" );
        
        }
        { //::SireMM::AmberParams::improperFunctions
        
            typedef ::SireMM::FourAtomFunctions ( ::SireMM::AmberParams::*improperFunctions_function_type)(  ) const;
            improperFunctions_function_type improperFunctions_function_value( &::SireMM::AmberParams::improperFunctions );
            
            AmberParams_exposer.def( 
                "improperFunctions"
                , improperFunctions_function_value
                , "Return all of the improper parameters converted to a set of FourAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::improperFunctions
        
            typedef ::SireMM::FourAtomFunctions ( ::SireMM::AmberParams::*improperFunctions_function_type)( ::SireCAS::Symbol const & ) const;
            improperFunctions_function_type improperFunctions_function_value( &::SireMM::AmberParams::improperFunctions );
            
            AmberParams_exposer.def( 
                "improperFunctions"
                , improperFunctions_function_value
                , ( bp::arg("PHI") )
                , "Return all of the improper parameters converted to a set of FourAtomFunctions" );
        
        }
        { //::SireMM::AmberParams::impropers
        
            typedef ::QHash< SireMol::ImproperID, QPair< SireMM::AmberDihedral, bool > > ( ::SireMM::AmberParams::*impropers_function_type)(  ) const;
            impropers_function_type impropers_function_value( &::SireMM::AmberParams::impropers );
            
            AmberParams_exposer.def( 
                "impropers"
                , impropers_function_value
                , "" );
        
        }
        { //::SireMM::AmberParams::info
        
            typedef ::SireMol::MoleculeInfo ( ::SireMM::AmberParams::*info_function_type)(  ) const;
            info_function_type info_function_value( &::SireMM::AmberParams::info );
            
            AmberParams_exposer.def( 
                "info"
                , info_function_value
                , "Return the layout of the molecule whose flexibility is contained\nin this object" );
        
        }
        { //::SireMM::AmberParams::isCompatibleWith
        
            typedef bool ( ::SireMM::AmberParams::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMM::AmberParams::isCompatibleWith );
            
            AmberParams_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , "Return whether or not this flexibility is compatible with the molecule\nwhose info is in molinfo" );
        
        }
        { //::SireMM::AmberParams::ljs
        
            typedef ::SireMM::AtomLJs ( ::SireMM::AmberParams::*ljs_function_type)(  ) const;
            ljs_function_type ljs_function_value( &::SireMM::AmberParams::ljs );
            
            AmberParams_exposer.def( 
                "ljs"
                , ljs_function_value
                , "Return the atom LJ parameters" );
        
        }
        { //::SireMM::AmberParams::masses
        
            typedef ::SireMol::AtomMasses ( ::SireMM::AmberParams::*masses_function_type)(  ) const;
            masses_function_type masses_function_value( &::SireMM::AmberParams::masses );
            
            AmberParams_exposer.def( 
                "masses"
                , masses_function_value
                , "Return the atom masses" );
        
        }
        { //::SireMM::AmberParams::nb14s
        
            typedef ::QHash< SireMol::BondID, SireMM::AmberNB14 > ( ::SireMM::AmberParams::*nb14s_function_type)(  ) const;
            nb14s_function_type nb14s_function_value( &::SireMM::AmberParams::nb14s );
            
            AmberParams_exposer.def( 
                "nb14s"
                , nb14s_function_value
                , "" );
        
        }
        AmberParams_exposer.def( bp::self != bp::self );
        AmberParams_exposer.def( bp::self + bp::self );
        { //::SireMM::AmberParams::operator=
        
            typedef ::SireMM::AmberParams & ( ::SireMM::AmberParams::*assign_function_type)( ::SireMM::AmberParams const & ) ;
            assign_function_type assign_function_value( &::SireMM::AmberParams::operator= );
            
            AmberParams_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AmberParams_exposer.def( bp::self == bp::self );
        { //::SireMM::AmberParams::propertyMap
        
            typedef ::SireBase::PropertyMap const & ( ::SireMM::AmberParams::*propertyMap_function_type)(  ) const;
            propertyMap_function_type propertyMap_function_value( &::SireMM::AmberParams::propertyMap );
            
            AmberParams_exposer.def( 
                "propertyMap"
                , propertyMap_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the property map that is used to find and update properties\nof the molecule" );
        
        }
        { //::SireMM::AmberParams::radiusSet
        
            typedef ::QString ( ::SireMM::AmberParams::*radiusSet_function_type)(  ) const;
            radiusSet_function_type radiusSet_function_value( &::SireMM::AmberParams::radiusSet );
            
            AmberParams_exposer.def( 
                "radiusSet"
                , radiusSet_function_value
                , "Return the radius set used by LEAP to assign the Born radii" );
        
        }
        { //::SireMM::AmberParams::remove
        
            typedef void ( ::SireMM::AmberParams::*remove_function_type)( ::SireMol::BondID const & ) ;
            remove_function_type remove_function_value( &::SireMM::AmberParams::remove );
            
            AmberParams_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("bond") )
                , "" );
        
        }
        { //::SireMM::AmberParams::remove
        
            typedef void ( ::SireMM::AmberParams::*remove_function_type)( ::SireMol::AngleID const & ) ;
            remove_function_type remove_function_value( &::SireMM::AmberParams::remove );
            
            AmberParams_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("angle") )
                , "" );
        
        }
        { //::SireMM::AmberParams::remove
        
            typedef void ( ::SireMM::AmberParams::*remove_function_type)( ::SireMol::DihedralID const & ) ;
            remove_function_type remove_function_value( &::SireMM::AmberParams::remove );
            
            AmberParams_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("dihedral") )
                , "" );
        
        }
        { //::SireMM::AmberParams::remove
        
            typedef void ( ::SireMM::AmberParams::*remove_function_type)( ::SireMol::ImproperID const & ) ;
            remove_function_type remove_function_value( &::SireMM::AmberParams::remove );
            
            AmberParams_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("improper") )
                , "" );
        
        }
        { //::SireMM::AmberParams::removeNB14
        
            typedef void ( ::SireMM::AmberParams::*removeNB14_function_type)( ::SireMol::BondID const & ) ;
            removeNB14_function_type removeNB14_function_value( &::SireMM::AmberParams::removeNB14 );
            
            AmberParams_exposer.def( 
                "removeNB14"
                , removeNB14_function_value
                , ( bp::arg("pair") )
                , "" );
        
        }
        { //::SireMM::AmberParams::setExcludedAtoms
        
            typedef void ( ::SireMM::AmberParams::*setExcludedAtoms_function_type)( ::SireMM::CLJNBPairs const & ) ;
            setExcludedAtoms_function_type setExcludedAtoms_function_value( &::SireMM::AmberParams::setExcludedAtoms );
            
            AmberParams_exposer.def( 
                "setExcludedAtoms"
                , setExcludedAtoms_function_value
                , ( bp::arg("excluded_atoms") )
                , "Set the excluded atoms of the molecule. This should be a\nCLJNBPairs with the value equal to 0 for atom0-atom1 pairs\nthat are excluded, and 1 for atom0-atom1 pairs that are\nto be included in the non-bonded calculation" );
        
        }
        { //::SireMM::AmberParams::setPropertyMap
        
            typedef void ( ::SireMM::AmberParams::*setPropertyMap_function_type)( ::SireBase::PropertyMap const & ) ;
            setPropertyMap_function_type setPropertyMap_function_value( &::SireMM::AmberParams::setPropertyMap );
            
            AmberParams_exposer.def( 
                "setPropertyMap"
                , setPropertyMap_function_value
                , ( bp::arg("map") )
                , "Set the property map that should be used to find and update properties\nof the molecule" );
        
        }
        { //::SireMM::AmberParams::setRadiusSet
        
            typedef void ( ::SireMM::AmberParams::*setRadiusSet_function_type)( ::QString const & ) ;
            setRadiusSet_function_type setRadiusSet_function_value( &::SireMM::AmberParams::setRadiusSet );
            
            AmberParams_exposer.def( 
                "setRadiusSet"
                , setRadiusSet_function_value
                , ( bp::arg("radius_set") )
                , "Set the radius set used by LEAP to assign the Born radii\nof the atoms. This is just a string that is used to label\nthe radius set in the PRM file" );
        
        }
        { //::SireMM::AmberParams::toString
        
            typedef ::QString ( ::SireMM::AmberParams::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::AmberParams::toString );
            
            AmberParams_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMM::AmberParams::treeChains
        
            typedef ::SireMol::AtomStringProperty ( ::SireMM::AmberParams::*treeChains_function_type)(  ) const;
            treeChains_function_type treeChains_function_value( &::SireMM::AmberParams::treeChains );
            
            AmberParams_exposer.def( 
                "treeChains"
                , treeChains_function_value
                , "Return all of the Amber treechain classification for all of the atoms" );
        
        }
        { //::SireMM::AmberParams::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::AmberParams::typeName );
            
            AmberParams_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::AmberParams::updateFrom
        
            typedef void ( ::SireMM::AmberParams::*updateFrom_function_type)( ::SireMol::MoleculeView const & ) ;
            updateFrom_function_type updateFrom_function_value( &::SireMM::AmberParams::updateFrom );
            
            AmberParams_exposer.def( 
                "updateFrom"
                , updateFrom_function_value
                , ( bp::arg("molview") )
                , "Update these parameters from the contents of the passed molecule. This\nwill only work if these parameters are compatible with this molecule" );
        
        }
        { //::SireMM::AmberParams::validate
        
            typedef ::QStringList ( ::SireMM::AmberParams::*validate_function_type)(  ) const;
            validate_function_type validate_function_value( &::SireMM::AmberParams::validate );
            
            AmberParams_exposer.def( 
                "validate"
                , validate_function_value
                , "Validate this set of parameters. This checks that all of the requirements\nfor an Amber set of parameters are met, e.g. that all Atom indicies are\ncontiguous and in-order, and that all atoms contiguously fill all residues\netc. This returns any errors as strings. An empty set of strings indicates\nthat there are no errors" );
        
        }
        { //::SireMM::AmberParams::validateAndFix
        
            typedef ::QStringList ( ::SireMM::AmberParams::*validateAndFix_function_type)(  ) ;
            validateAndFix_function_type validateAndFix_function_value( &::SireMM::AmberParams::validateAndFix );
            
            AmberParams_exposer.def( 
                "validateAndFix"
                , validateAndFix_function_value
                , "" );
        
        }
        AmberParams_exposer.staticmethod( "typeName" );
        AmberParams_exposer.def( "__copy__", &__copy__);
        AmberParams_exposer.def( "__deepcopy__", &__copy__);
        AmberParams_exposer.def( "clone", &__copy__);
        AmberParams_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::AmberParams >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AmberParams_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::AmberParams >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AmberParams_exposer.def( "__str__", &__str__< ::SireMM::AmberParams > );
        AmberParams_exposer.def( "__repr__", &__str__< ::SireMM::AmberParams > );
    }

}
