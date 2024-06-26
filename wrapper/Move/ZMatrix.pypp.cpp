// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ZMatrix.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMol/angleid.h"

#include "SireMol/atommatcher.h"

#include "SireMol/bondid.h"

#include "SireMol/dihedralid.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireMove/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/convert.h"

#include "SireUnits/units.h"

#include "zmatrix.h"

#include <QDebug>

#include <QElapsedTimer>

#include <QTime>

#include "zmatrix.h"

SireMove::ZMatrix __copy__(const SireMove::ZMatrix &other){ return SireMove::ZMatrix(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_ZMatrix_class(){

    { //::SireMove::ZMatrix
        typedef bp::class_< SireMove::ZMatrix, bp::bases< SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > ZMatrix_exposer_t;
        ZMatrix_exposer_t ZMatrix_exposer = ZMatrix_exposer_t( "ZMatrix", "This class holds the z-matrix of a molecule. A z-matrix\nis be used to provide an internal coordinates representation\nof a molecule, so that internal moves may be performed\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope ZMatrix_scope( ZMatrix_exposer );
        ZMatrix_exposer.def( bp::init< SireMol::Molecule const & >(( bp::arg("molecule") ), "Construct to hold the z-matrix for the passed molecule") );
        ZMatrix_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "Construct to hold the z-matrix for the molecule whose\nlayout information is held in molinfo") );
        ZMatrix_exposer.def( bp::init< SireMove::ZMatrix const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::ZMatrix::add
        
            typedef void ( ::SireMove::ZMatrix::*add_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) ;
            add_function_type add_function_value( &::SireMove::ZMatrix::add );
            
            ZMatrix_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle"), bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Add a line to the z-matrix that gives the coordinates of the\natom atom based on the passed bond, angle and dihedal atoms\nAn exception is raised if adding this line would lead\nto a circular reference\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::add
        
            typedef void ( ::SireMove::ZMatrix::*add_function_type)( ::SireMol::DihedralID const & ) ;
            add_function_type add_function_value( &::SireMove::ZMatrix::add );
            
            ZMatrix_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Add a line to the z-matrix that gives the coordinates of the\npassed dihedral\nAn exception is raised if adding this line would lead\nto a circular reference\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::add
        
            typedef void ( ::SireMove::ZMatrix::*add_function_type)( ::SireMove::ZMatrixLine const & ) ;
            add_function_type add_function_value( &::SireMove::ZMatrix::add );
            
            ZMatrix_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("zmatline") )
                , bp::release_gil_policy()
                , "Add the z-matrix line zmatline.\nAn exception is raised if adding this line would lead to\na circular reference\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::angleDelta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrix::*angleDelta_function_type)( ::SireMol::AtomID const & ) const;
            angleDelta_function_type angleDelta_function_value( &::SireMove::ZMatrix::angleDelta );
            
            ZMatrix_exposer.def( 
                "angleDelta"
                , angleDelta_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Return the maximum amount that the angle to atom atom\nshould be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::angleDelta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrix::*angleDelta_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            angleDelta_function_type angleDelta_function_value( &::SireMove::ZMatrix::angleDelta );
            
            ZMatrix_exposer.def( 
                "angleDelta"
                , angleDelta_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle") )
                , bp::release_gil_policy()
                , "Return the maximum amount that the angle between atoms\natom-bond-angle should be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::at
        
            typedef ::SireMove::ZMatrixLine const & ( ::SireMove::ZMatrix::*at_function_type)( ::SireMol::AtomID const & ) const;
            at_function_type at_function_value( &::SireMove::ZMatrix::at );
            
            ZMatrix_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("atom") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the line for the atom identified by atom. This\nraises an exception if this is not a valid ID, or if\nthere is no z-matrix line for this atom\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::atomBuildOrder
        
            typedef ::QVector< int > const & ( ::SireMove::ZMatrix::*atomBuildOrder_function_type)(  ) const;
            atomBuildOrder_function_type atomBuildOrder_function_value( &::SireMove::ZMatrix::atomBuildOrder );
            
            ZMatrix_exposer.def( 
                "atomBuildOrder"
                , atomBuildOrder_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the correct atom build order for this z-matrix" );
        
        }
        { //::SireMove::ZMatrix::bondDelta
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::ZMatrix::*bondDelta_function_type)( ::SireMol::AtomID const & ) const;
            bondDelta_function_type bondDelta_function_value( &::SireMove::ZMatrix::bondDelta );
            
            ZMatrix_exposer.def( 
                "bondDelta"
                , bondDelta_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Return the maximum amount that the bond to atom atom\nshould be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::bondDelta
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::ZMatrix::*bondDelta_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            bondDelta_function_type bondDelta_function_value( &::SireMove::ZMatrix::bondDelta );
            
            ZMatrix_exposer.def( 
                "bondDelta"
                , bondDelta_function_value
                , ( bp::arg("atom"), bp::arg("bond") )
                , bp::release_gil_policy()
                , "Return the maximum amount that the bond between atoms\natom-bond should be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::contains
        
            typedef bool ( ::SireMove::ZMatrix::*contains_function_type)( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMove::ZMatrix::contains );
            
            ZMatrix_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix contains a line for the\natom with ID atom\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::contains
        
            typedef bool ( ::SireMove::ZMatrix::*contains_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMove::ZMatrix::contains );
            
            ZMatrix_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atom"), bp::arg("bond") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix contains the bond between\natoms atom-bond\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::contains
        
            typedef bool ( ::SireMove::ZMatrix::*contains_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMove::ZMatrix::contains );
            
            ZMatrix_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix contains the angle\natom-bond-angle\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::contains
        
            typedef bool ( ::SireMove::ZMatrix::*contains_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMove::ZMatrix::contains );
            
            ZMatrix_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle"), bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix contains the dihedral\natom-bond-angle-dihedral\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::contains
        
            typedef bool ( ::SireMove::ZMatrix::*contains_function_type)( ::SireMol::BondID const & ) const;
            contains_function_type contains_function_value( &::SireMove::ZMatrix::contains );
            
            ZMatrix_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("bond") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix contains bond\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::contains
        
            typedef bool ( ::SireMove::ZMatrix::*contains_function_type)( ::SireMol::AngleID const & ) const;
            contains_function_type contains_function_value( &::SireMove::ZMatrix::contains );
            
            ZMatrix_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("angle") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix contains angle\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::contains
        
            typedef bool ( ::SireMove::ZMatrix::*contains_function_type)( ::SireMol::DihedralID const & ) const;
            contains_function_type contains_function_value( &::SireMove::ZMatrix::contains );
            
            ZMatrix_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix contains dihedral\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::delta
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::ZMatrix::*delta_function_type)( ::SireMol::BondID const & ) const;
            delta_function_type delta_function_value( &::SireMove::ZMatrix::delta );
            
            ZMatrix_exposer.def( 
                "delta"
                , delta_function_value
                , ( bp::arg("bond") )
                , bp::release_gil_policy()
                , "Return the maximum amount that bond should be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::delta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrix::*delta_function_type)( ::SireMol::AngleID const & ) const;
            delta_function_type delta_function_value( &::SireMove::ZMatrix::delta );
            
            ZMatrix_exposer.def( 
                "delta"
                , delta_function_value
                , ( bp::arg("angle") )
                , bp::release_gil_policy()
                , "Return the maximum amount that angle should be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::delta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrix::*delta_function_type)( ::SireMol::DihedralID const & ) const;
            delta_function_type delta_function_value( &::SireMove::ZMatrix::delta );
            
            ZMatrix_exposer.def( 
                "delta"
                , delta_function_value
                , ( bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Return the maximum amount that dihedral should be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::dihedralDelta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrix::*dihedralDelta_function_type)( ::SireMol::AtomID const & ) const;
            dihedralDelta_function_type dihedralDelta_function_value( &::SireMove::ZMatrix::dihedralDelta );
            
            ZMatrix_exposer.def( 
                "dihedralDelta"
                , dihedralDelta_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Return the maximum amount that the dihedral to atom atom\nshould be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::dihedralDelta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::ZMatrix::*dihedralDelta_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            dihedralDelta_function_type dihedralDelta_function_value( &::SireMove::ZMatrix::dihedralDelta );
            
            ZMatrix_exposer.def( 
                "dihedralDelta"
                , dihedralDelta_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle"), bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Return the maximum amount that the dihedral between atoms\natom-bond-angle-dihedral should be changed\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::getIndex
        
            typedef int ( ::SireMove::ZMatrix::*getIndex_function_type)( ::SireMol::AtomID const & ) const;
            getIndex_function_type getIndex_function_value( &::SireMove::ZMatrix::getIndex );
            
            ZMatrix_exposer.def( 
                "getIndex"
                , getIndex_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Return the index of the z-matrix line that positions the atom\nwith ID atom\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::getIndex
        
            typedef int ( ::SireMove::ZMatrix::*getIndex_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            getIndex_function_type getIndex_function_value( &::SireMove::ZMatrix::getIndex );
            
            ZMatrix_exposer.def( 
                "getIndex"
                , getIndex_function_value
                , ( bp::arg("atom"), bp::arg("bond") )
                , bp::release_gil_policy()
                , "Return the index of the z-matrix line that positions the\nbond between atoms atom-bond\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::getIndex
        
            typedef int ( ::SireMove::ZMatrix::*getIndex_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            getIndex_function_type getIndex_function_value( &::SireMove::ZMatrix::getIndex );
            
            ZMatrix_exposer.def( 
                "getIndex"
                , getIndex_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle") )
                , bp::release_gil_policy()
                , "Return the index of the z-matrix line that positions the\nbond between atoms atom-bond-angle\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::getIndex
        
            typedef int ( ::SireMove::ZMatrix::*getIndex_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            getIndex_function_type getIndex_function_value( &::SireMove::ZMatrix::getIndex );
            
            ZMatrix_exposer.def( 
                "getIndex"
                , getIndex_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle"), bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Return the index of the z-matrix line that positions the\nbond between atoms atom-bond-angle-dihedral\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::getIndex
        
            typedef int ( ::SireMove::ZMatrix::*getIndex_function_type)( ::SireMol::BondID const & ) const;
            getIndex_function_type getIndex_function_value( &::SireMove::ZMatrix::getIndex );
            
            ZMatrix_exposer.def( 
                "getIndex"
                , getIndex_function_value
                , ( bp::arg("bond") )
                , bp::release_gil_policy()
                , "Return the index of the z-matrix line that defines the\nbond bond\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::getIndex
        
            typedef int ( ::SireMove::ZMatrix::*getIndex_function_type)( ::SireMol::AngleID const & ) const;
            getIndex_function_type getIndex_function_value( &::SireMove::ZMatrix::getIndex );
            
            ZMatrix_exposer.def( 
                "getIndex"
                , getIndex_function_value
                , ( bp::arg("angle") )
                , bp::release_gil_policy()
                , "Return the index of the z-matrix line that defines the\nangle angle\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::getIndex
        
            typedef int ( ::SireMove::ZMatrix::*getIndex_function_type)( ::SireMol::DihedralID const & ) const;
            getIndex_function_type getIndex_function_value( &::SireMove::ZMatrix::getIndex );
            
            ZMatrix_exposer.def( 
                "getIndex"
                , getIndex_function_value
                , ( bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Return the index of the z-matrix line that defines the\ndihedral dihedral\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::index
        
            typedef ::QHash< SireMol::AtomIdx, int > const & ( ::SireMove::ZMatrix::*index_function_type)(  ) const;
            index_function_type index_function_value( &::SireMove::ZMatrix::index );
            
            ZMatrix_exposer.def( 
                "index"
                , index_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the index of AtomIdx to z-matrix line number. This\nis used to index the output of ZMatrix::lines()" );
        
        }
        { //::SireMove::ZMatrix::info
        
            typedef ::SireMol::MoleculeInfoData const & ( ::SireMove::ZMatrix::*info_function_type)(  ) const;
            info_function_type info_function_value( &::SireMove::ZMatrix::info );
            
            ZMatrix_exposer.def( 
                "info"
                , info_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the layout of the molecule whose z-matrix is contained\nin this object" );
        
        }
        { //::SireMove::ZMatrix::isCompatibleWith
        
            typedef bool ( ::SireMove::ZMatrix::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMove::ZMatrix::isCompatibleWith );
            
            ZMatrix_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "Return whether or not this z-matrix is compatible with the\nthe molecule whose info is in molinfo" );
        
        }
        { //::SireMove::ZMatrix::lines
        
            typedef ::QVector< SireMove::ZMatrixLine > const & ( ::SireMove::ZMatrix::*lines_function_type)(  ) const;
            lines_function_type lines_function_value( &::SireMove::ZMatrix::lines );
            
            ZMatrix_exposer.def( 
                "lines"
                , lines_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return all of the lines of the z-matrix" );
        
        }
        { //::SireMove::ZMatrix::matchToSelection
        
            typedef ::SireMove::ZMatrix ( ::SireMove::ZMatrix::*matchToSelection_function_type)( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMove::ZMatrix::matchToSelection );
            
            ZMatrix_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "Return a z-matrix that only contains lines that involve the atoms\nthat are in selection" );
        
        }
        { //::SireMove::ZMatrix::merge
        
            typedef ::SireBase::PropertyList ( ::SireMove::ZMatrix::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMove::ZMatrix::merge );
            
            ZMatrix_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMove::ZMatrix::nLines
        
            typedef int ( ::SireMove::ZMatrix::*nLines_function_type)(  ) const;
            nLines_function_type nLines_function_value( &::SireMove::ZMatrix::nLines );
            
            ZMatrix_exposer.def( 
                "nLines"
                , nLines_function_value
                , bp::release_gil_policy()
                , "Return the number of lines in this z-matrix" );
        
        }
        ZMatrix_exposer.def( bp::self != bp::self );
        { //::SireMove::ZMatrix::operator=
        
            typedef ::SireMove::ZMatrix & ( ::SireMove::ZMatrix::*assign_function_type)( ::SireMove::ZMatrix const & ) ;
            assign_function_type assign_function_value( &::SireMove::ZMatrix::operator= );
            
            ZMatrix_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ZMatrix_exposer.def( bp::self == bp::self );
        { //::SireMove::ZMatrix::operator[]
        
            typedef ::SireMove::ZMatrixLine const & ( ::SireMove::ZMatrix::*__getitem___function_type)( ::SireMol::AtomID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMove::ZMatrix::operator[] );
            
            ZMatrix_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("atom") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMove::ZMatrix::remove
        
            typedef void ( ::SireMove::ZMatrix::*remove_function_type)( ::SireMol::AtomID const & ) ;
            remove_function_type remove_function_value( &::SireMove::ZMatrix::remove );
            
            ZMatrix_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Remove the z-matrix line that gives the coordinates of the\natom atom. This does nothing if there is no z-matrix\nline for this atom - note this removes the lines for\nall matching atoms\nThrow: SireMol::missing_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::remove
        
            typedef void ( ::SireMove::ZMatrix::*remove_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) ;
            remove_function_type remove_function_value( &::SireMove::ZMatrix::remove );
            
            ZMatrix_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle"), bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Remove the z-matrix line involving the specified atom, bond,\nangle and dihedral. This matches only a single line, and\nwill do nothing if this line is not in this z-matrix\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::remove
        
            typedef void ( ::SireMove::ZMatrix::*remove_function_type)( ::SireMol::DihedralID const & ) ;
            remove_function_type remove_function_value( &::SireMove::ZMatrix::remove );
            
            ZMatrix_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("dihedral") )
                , bp::release_gil_policy()
                , "Remove the z-matrix line involving the specified dihedral.\nThis matches only a single line, and\nwill do nothing if this line is not in this z-matrix\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::ZMatrix::remove
        
            typedef void ( ::SireMove::ZMatrix::*remove_function_type)( ::SireMove::ZMatrixLine const & ) ;
            remove_function_type remove_function_value( &::SireMove::ZMatrix::remove );
            
            ZMatrix_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("zmatline") )
                , bp::release_gil_policy()
                , "Remove the z-matrix line zmatline from this z-matrix" );
        
        }
        { //::SireMove::ZMatrix::setAngleDelta
        
            typedef void ( ::SireMove::ZMatrix::*setAngleDelta_function_type)( ::SireMol::AtomID const &,::SireUnits::Dimension::Angle const & ) ;
            setAngleDelta_function_type setAngleDelta_function_value( &::SireMove::ZMatrix::setAngleDelta );
            
            ZMatrix_exposer.def( 
                "setAngleDelta"
                , setAngleDelta_function_value
                , ( bp::arg("atom"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the angle for the atom atom\ncan be changed to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setAngleDelta
        
            typedef void ( ::SireMove::ZMatrix::*setAngleDelta_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireUnits::Dimension::Angle const & ) ;
            setAngleDelta_function_type setAngleDelta_function_value( &::SireMove::ZMatrix::setAngleDelta );
            
            ZMatrix_exposer.def( 
                "setAngleDelta"
                , setAngleDelta_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the angle between atoms atom-bond-angle\ncan be changed by to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setBondDelta
        
            typedef void ( ::SireMove::ZMatrix::*setBondDelta_function_type)( ::SireMol::AtomID const &,::SireUnits::Dimension::Length const & ) ;
            setBondDelta_function_type setBondDelta_function_value( &::SireMove::ZMatrix::setBondDelta );
            
            ZMatrix_exposer.def( 
                "setBondDelta"
                , setBondDelta_function_value
                , ( bp::arg("atom"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the bond for the atom atom\ncan be moved to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setBondDelta
        
            typedef void ( ::SireMove::ZMatrix::*setBondDelta_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireUnits::Dimension::Length const & ) ;
            setBondDelta_function_type setBondDelta_function_value( &::SireMove::ZMatrix::setBondDelta );
            
            ZMatrix_exposer.def( 
                "setBondDelta"
                , setBondDelta_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the bond between atoms atom-bond\ncan be changed by to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setDelta
        
            typedef void ( ::SireMove::ZMatrix::*setDelta_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Length const & ) ;
            setDelta_function_type setDelta_function_value( &::SireMove::ZMatrix::setDelta );
            
            ZMatrix_exposer.def( 
                "setDelta"
                , setDelta_function_value
                , ( bp::arg("bond"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the bond bond\ncan be changed by to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setDelta
        
            typedef void ( ::SireMove::ZMatrix::*setDelta_function_type)( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle const & ) ;
            setDelta_function_type setDelta_function_value( &::SireMove::ZMatrix::setDelta );
            
            ZMatrix_exposer.def( 
                "setDelta"
                , setDelta_function_value
                , ( bp::arg("angle"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the angle angle\ncan be changed by to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setDelta
        
            typedef void ( ::SireMove::ZMatrix::*setDelta_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle const & ) ;
            setDelta_function_type setDelta_function_value( &::SireMove::ZMatrix::setDelta );
            
            ZMatrix_exposer.def( 
                "setDelta"
                , setDelta_function_value
                , ( bp::arg("dihedral"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the dihedral dihedral\ncan be changed by to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setDihedralDelta
        
            typedef void ( ::SireMove::ZMatrix::*setDihedralDelta_function_type)( ::SireMol::AtomID const &,::SireUnits::Dimension::Angle const & ) ;
            setDihedralDelta_function_type setDihedralDelta_function_value( &::SireMove::ZMatrix::setDihedralDelta );
            
            ZMatrix_exposer.def( 
                "setDihedralDelta"
                , setDihedralDelta_function_value
                , ( bp::arg("atom"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the dihedral for the atom atom\ncan be changed to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::setDihedralDelta
        
            typedef void ( ::SireMove::ZMatrix::*setDihedralDelta_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireUnits::Dimension::Angle const & ) ;
            setDihedralDelta_function_type setDihedralDelta_function_value( &::SireMove::ZMatrix::setDihedralDelta );
            
            ZMatrix_exposer.def( 
                "setDihedralDelta"
                , setDihedralDelta_function_value
                , ( bp::arg("atom"), bp::arg("bond"), bp::arg("angle"), bp::arg("dihedral"), bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum amount that the dihedral between atoms\natom-bond-angle-dihedral can be changed by to delta\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\nThrow: SireMove::zmatrix_error\n" );
        
        }
        { //::SireMove::ZMatrix::toString
        
            typedef ::QString ( ::SireMove::ZMatrix::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::ZMatrix::toString );
            
            ZMatrix_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this z-matrix" );
        
        }
        { //::SireMove::ZMatrix::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::ZMatrix::typeName );
            
            ZMatrix_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ZMatrix_exposer.staticmethod( "typeName" );
        ZMatrix_exposer.def( "__copy__", &__copy__);
        ZMatrix_exposer.def( "__deepcopy__", &__copy__);
        ZMatrix_exposer.def( "clone", &__copy__);
        ZMatrix_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::ZMatrix >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ZMatrix_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::ZMatrix >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ZMatrix_exposer.def_pickle(sire_pickle_suite< ::SireMove::ZMatrix >());
        ZMatrix_exposer.def( "__str__", &__str__< ::SireMove::ZMatrix > );
        ZMatrix_exposer.def( "__repr__", &__str__< ::SireMove::ZMatrix > );
    }

}
