// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "TriclinicBox.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/align.h"

#include "SireMaths/rangenerator.h"

#include "SireStream/datastream.h"

#include "coordgroup.h"

#include "triclinicbox.h"

#include <QDebug>

#include <QList>

#include <QPair>

#include <cmath>

#include <limits>

#include "triclinicbox.h"

SireVol::TriclinicBox __copy__(const SireVol::TriclinicBox &other){ return SireVol::TriclinicBox(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_TriclinicBox_class(){

    { //::SireVol::TriclinicBox
        typedef bp::class_< SireVol::TriclinicBox, bp::bases< SireVol::Cartesian, SireVol::Space, SireBase::Property > > TriclinicBox_exposer_t;
        TriclinicBox_exposer_t TriclinicBox_exposer = TriclinicBox_exposer_t( "TriclinicBox", "\nA TriclinicBox is a volume that represents standard periodic boundary conditions\n(a 3D box replicated to infinity along all three dimensions).\n\nTo support triclinic boxes that work across a range of molecular simulation\nengines, e.g. AMBER, GROMACS, OpenMM, we represent the triclinic space in\nreduced form, using the approach documented in Appendix A of Chapter 3 from\nMolecular dynamics of sense and sensibility in processing and analysis of data\nby Tsjerk A. Wassenaar.\n\nAuthor: Lester Hedges\n", bp::init< >("Construct a default TriclinicBox (large volume)") );
        bp::scope TriclinicBox_scope( TriclinicBox_exposer );
        TriclinicBox_exposer.def( bp::init< SireMaths::Vector const &, SireMaths::Vector const &, SireMaths::Vector const &, bp::optional< bool, bool > >(( bp::arg("v0"), bp::arg("v1"), bp::arg("v2"), bp::arg("auto_rotate")=(bool)(false), bp::arg("auto_reduce")=(bool)(false) ), "Construct a triclinic box from lattice vectors.\n\nPar:am v0\nThe first lattice vector.\n\nPar:am v1\nThe second lattice vector.\n\nPar:am v2\nThe third lattice vector.\n\nPar:am auto_rotate\nWhether to automatically rotate the box to comply with the\nconstraints of molecular dynamics engines, i.e. vector0 aligned\nwith x axis, vector1 in x-y plane, and vector2 with positive\nz component.\n\nPar:am auto_reduce\nWhether to automatically perform a lattice reduction on the\nbox.\n") );
        TriclinicBox_exposer.def( bp::init< double, double, double, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, bp::optional< bool, bool > >(( bp::arg("a"), bp::arg("b"), bp::arg("c"), bp::arg("alpha"), bp::arg("beta"), bp::arg("gamma"), bp::arg("auto_rotate")=(bool)(false), bp::arg("auto_reduce")=(bool)(false) ), "Construct a triclinic box from box lengths and angles.\n\nPar:am a\nThe length of the first box vector.\n\nPar:am b\nThe length of the second box vector.\n\nPar:am c\nThe length of the third box vector.\n\nPar:am alpha\nThe angle between the second and third box vectors.\n\nPar:am beta\nThe angle between the first and third box vectors.\n\nPar:am gamma\nThe angle between the second and first box vectors.\n\nPar:am auto_rotate\nWhether to automatically rotate the box to comply with the\nconstraints of molecular dynamics engines, i.e. vector0 aligned\nwith x axis, vector1 in x-y plane, and vector2 with positive\nz component:w.\n\nPar:am auto_reduce\nWhether to automatically perform a lattice reduction on the\nbox.\n") );
        TriclinicBox_exposer.def( bp::init< SireVol::TriclinicBox const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireVol::TriclinicBox::alpha
        
            typedef double ( ::SireVol::TriclinicBox::*alpha_function_type)(  ) const;
            alpha_function_type alpha_function_value( &::SireVol::TriclinicBox::alpha );
            
            TriclinicBox_exposer.def( 
                "alpha"
                , alpha_function_value
                , bp::release_gil_policy()
                , "Return the angle between v1 and v2 in degrees." );
        
        }
        { //::SireVol::TriclinicBox::beta
        
            typedef double ( ::SireVol::TriclinicBox::*beta_function_type)(  ) const;
            beta_function_type beta_function_value( &::SireVol::TriclinicBox::beta );
            
            TriclinicBox_exposer.def( 
                "beta"
                , beta_function_value
                , bp::release_gil_policy()
                , "Return the angle between v0 and v2 in degrees." );
        
        }
        { //::SireVol::TriclinicBox::beyond
        
            typedef bool ( ::SireVol::TriclinicBox::*beyond_function_type)( double,::SireVol::AABox const &,::SireVol::AABox const & ) const;
            beyond_function_type beyond_function_value( &::SireVol::TriclinicBox::beyond );
            
            TriclinicBox_exposer.def( 
                "beyond"
                , beyond_function_value
                , ( bp::arg("dist"), bp::arg("aabox0"), bp::arg("aabox1") )
                , bp::release_gil_policy()
                , "Return whether or not two groups enclosed by the AABoxes aabox0 and\naabox1 are definitely beyond the cutoff distance dist" );
        
        }
        { //::SireVol::TriclinicBox::beyond
        
            typedef bool ( ::SireVol::TriclinicBox::*beyond_function_type)( double,::SireVol::CoordGroup const &,::SireVol::CoordGroup const & ) const;
            beyond_function_type beyond_function_value( &::SireVol::TriclinicBox::beyond );
            
            TriclinicBox_exposer.def( 
                "beyond"
                , beyond_function_value
                , ( bp::arg("dist"), bp::arg("group0"), bp::arg("group1") )
                , bp::release_gil_policy()
                , "Return whether or not these two groups are definitely beyond the cutoff distance." );
        
        }
        { //::SireVol::TriclinicBox::boxMatrix
        
            typedef ::SireMaths::Matrix ( ::SireVol::TriclinicBox::*boxMatrix_function_type)(  ) const;
            boxMatrix_function_type boxMatrix_function_value( &::SireVol::TriclinicBox::boxMatrix );
            
            TriclinicBox_exposer.def( 
                "boxMatrix"
                , boxMatrix_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::TriclinicBox::calcAngle
        
            typedef ::SireUnits::Dimension::Angle ( ::SireVol::TriclinicBox::*calcAngle_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const &,::SireMaths::Vector const & ) const;
            calcAngle_function_type calcAngle_function_value( &::SireVol::TriclinicBox::calcAngle );
            
            TriclinicBox_exposer.def( 
                "calcAngle"
                , calcAngle_function_value
                , ( bp::arg("point0"), bp::arg("point1"), bp::arg("point2") )
                , bp::release_gil_policy()
                , "Calculate the angle between the passed three points. This should return\nthe acute angle between the points, which should lie between 0 and 180 degrees" );
        
        }
        { //::SireVol::TriclinicBox::calcDihedral
        
            typedef ::SireUnits::Dimension::Angle ( ::SireVol::TriclinicBox::*calcDihedral_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const &,::SireMaths::Vector const &,::SireMaths::Vector const & ) const;
            calcDihedral_function_type calcDihedral_function_value( &::SireVol::TriclinicBox::calcDihedral );
            
            TriclinicBox_exposer.def( 
                "calcDihedral"
                , calcDihedral_function_value
                , ( bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("point3") )
                , bp::release_gil_policy()
                , "Calculate the torsion angle between the passed four points. This should\nreturn the torsion angle measured clockwise when looking down the\ntorsion from point0-point1-point2-point3. This will lie between 0 and 360\ndegrees" );
        
        }
        { //::SireVol::TriclinicBox::calcDist
        
            typedef double ( ::SireVol::TriclinicBox::*calcDist_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const & ) const;
            calcDist_function_type calcDist_function_value( &::SireVol::TriclinicBox::calcDist );
            
            TriclinicBox_exposer.def( 
                "calcDist"
                , calcDist_function_value
                , ( bp::arg("point0"), bp::arg("point1") )
                , bp::release_gil_policy()
                , "Calculate the distance between two points" );
        
        }
        { //::SireVol::TriclinicBox::calcDist
        
            typedef double ( ::SireVol::TriclinicBox::*calcDist_function_type)( ::SireVol::CoordGroup const &,::SireVol::CoordGroup const &,::SireVol::DistMatrix & ) const;
            calcDist_function_type calcDist_function_value( &::SireVol::TriclinicBox::calcDist );
            
            TriclinicBox_exposer.def( 
                "calcDist"
                , calcDist_function_value
                , ( bp::arg("group1"), bp::arg("group2"), bp::arg("distmat") )
                , bp::release_gil_policy()
                , "Populate the matrix mat with the distances between all of the\natoms of the two CoordGroups. Return the shortest distance^2 between the two\nCoordGroups." );
        
        }
        { //::SireVol::TriclinicBox::calcDist
        
            typedef double ( ::SireVol::TriclinicBox::*calcDist_function_type)( ::SireVol::CoordGroup const &,::SireMaths::Vector const &,::SireVol::DistMatrix & ) const;
            calcDist_function_type calcDist_function_value( &::SireVol::TriclinicBox::calcDist );
            
            TriclinicBox_exposer.def( 
                "calcDist"
                , calcDist_function_value
                , ( bp::arg("group"), bp::arg("point"), bp::arg("mat") )
                , bp::release_gil_policy()
                , "Populate the matrix mat with the distances between all of the\natoms of the passed CoordGroup to the passed point. Return the shortest\ndistance." );
        
        }
        { //::SireVol::TriclinicBox::calcDist2
        
            typedef double ( ::SireVol::TriclinicBox::*calcDist2_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const & ) const;
            calcDist2_function_type calcDist2_function_value( &::SireVol::TriclinicBox::calcDist2 );
            
            TriclinicBox_exposer.def( 
                "calcDist2"
                , calcDist2_function_value
                , ( bp::arg("point0"), bp::arg("point1") )
                , bp::release_gil_policy()
                , "Calculate the distance squared between two points" );
        
        }
        { //::SireVol::TriclinicBox::calcDist2
        
            typedef double ( ::SireVol::TriclinicBox::*calcDist2_function_type)( ::SireVol::CoordGroup const &,::SireMaths::Vector const &,::SireVol::DistMatrix & ) const;
            calcDist2_function_type calcDist2_function_value( &::SireVol::TriclinicBox::calcDist2 );
            
            TriclinicBox_exposer.def( 
                "calcDist2"
                , calcDist2_function_value
                , ( bp::arg("group"), bp::arg("point"), bp::arg("mat") )
                , bp::release_gil_policy()
                , "Populate the matrix mat with the distances squared between all of the\natoms of the passed CoordGroup to the passed point. Return the shortest\ndistance." );
        
        }
        { //::SireVol::TriclinicBox::calcDist2
        
            typedef double ( ::SireVol::TriclinicBox::*calcDist2_function_type)( ::SireVol::CoordGroup const &,::SireVol::CoordGroup const &,::SireVol::DistMatrix & ) const;
            calcDist2_function_type calcDist2_function_value( &::SireVol::TriclinicBox::calcDist2 );
            
            TriclinicBox_exposer.def( 
                "calcDist2"
                , calcDist2_function_value
                , ( bp::arg("group1"), bp::arg("group2"), bp::arg("distmat") )
                , bp::release_gil_policy()
                , "Populate the matrix mat with the distances^2 between all of the\natoms of the two CoordGroups. Return the shortest distance between the\ntwo CoordGroups." );
        
        }
        { //::SireVol::TriclinicBox::calcDistVector
        
            typedef ::SireMaths::DistVector ( ::SireVol::TriclinicBox::*calcDistVector_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const & ) const;
            calcDistVector_function_type calcDistVector_function_value( &::SireVol::TriclinicBox::calcDistVector );
            
            TriclinicBox_exposer.def( 
                "calcDistVector"
                , calcDistVector_function_value
                , ( bp::arg("point0"), bp::arg("point1") )
                , bp::release_gil_policy()
                , "Calculate the distance vector between two points" );
        
        }
        { //::SireVol::TriclinicBox::calcDistVectors
        
            typedef double ( ::SireVol::TriclinicBox::*calcDistVectors_function_type)( ::SireVol::CoordGroup const &,::SireVol::CoordGroup const &,::SireVol::DistVectorMatrix & ) const;
            calcDistVectors_function_type calcDistVectors_function_value( &::SireVol::TriclinicBox::calcDistVectors );
            
            TriclinicBox_exposer.def( 
                "calcDistVectors"
                , calcDistVectors_function_value
                , ( bp::arg("group1"), bp::arg("group2"), bp::arg("distmat") )
                , bp::release_gil_policy()
                , "Populate the matrix distmat between all the points of the two CoordGroups\ngroup1 and group2 - the returned matrix has the vectors pointing\nfrom each point in group1 to each point in group2. This returns\nthe shortest distance between two points in the group" );
        
        }
        { //::SireVol::TriclinicBox::calcDistVectors
        
            typedef double ( ::SireVol::TriclinicBox::*calcDistVectors_function_type)( ::SireVol::CoordGroup const &,::SireMaths::Vector const &,::SireVol::DistVectorMatrix & ) const;
            calcDistVectors_function_type calcDistVectors_function_value( &::SireVol::TriclinicBox::calcDistVectors );
            
            TriclinicBox_exposer.def( 
                "calcDistVectors"
                , calcDistVectors_function_value
                , ( bp::arg("group"), bp::arg("point"), bp::arg("distmat") )
                , bp::release_gil_policy()
                , "Populate the matrix distmat between all the points passed CoordGroup\nto the point point - the returned matrix has the vectors pointing\nfrom the point to each point in group. This returns\nthe shortest distance." );
        
        }
        { //::SireVol::TriclinicBox::calcInvDist
        
            typedef double ( ::SireVol::TriclinicBox::*calcInvDist_function_type)( ::SireVol::CoordGroup const &,::SireVol::CoordGroup const &,::SireVol::DistMatrix & ) const;
            calcInvDist_function_type calcInvDist_function_value( &::SireVol::TriclinicBox::calcInvDist );
            
            TriclinicBox_exposer.def( 
                "calcInvDist"
                , calcInvDist_function_value
                , ( bp::arg("group1"), bp::arg("group2"), bp::arg("distmat") )
                , bp::release_gil_policy()
                , "Populate the matrix mat with the inverse distances between all of the\natoms of the two CoordGroups. Return the shortest distance between the two CoordGroups." );
        
        }
        { //::SireVol::TriclinicBox::calcInvDist2
        
            typedef double ( ::SireVol::TriclinicBox::*calcInvDist2_function_type)( ::SireVol::CoordGroup const &,::SireVol::CoordGroup const &,::SireVol::DistMatrix & ) const;
            calcInvDist2_function_type calcInvDist2_function_value( &::SireVol::TriclinicBox::calcInvDist2 );
            
            TriclinicBox_exposer.def( 
                "calcInvDist2"
                , calcInvDist2_function_value
                , ( bp::arg("group1"), bp::arg("group2"), bp::arg("distmat") )
                , bp::release_gil_policy()
                , "Populate the matrix mat with the inverse distances^2 between all of the\natoms of the two CoordGroups. Return the shortest distance between the two CoordGroups." );
        
        }
        { //::SireVol::TriclinicBox::cellMatrix
        
            typedef ::SireMaths::Matrix ( ::SireVol::TriclinicBox::*cellMatrix_function_type)(  ) const;
            cellMatrix_function_type cellMatrix_function_value( &::SireVol::TriclinicBox::cellMatrix );
            
            TriclinicBox_exposer.def( 
                "cellMatrix"
                , cellMatrix_function_value
                , bp::release_gil_policy()
                , "Return the cell matrix." );
        
        }
        { //::SireVol::TriclinicBox::cubic
        
            typedef ::SireVol::TriclinicBox ( *cubic_function_type )( double );
            cubic_function_type cubic_function_value( &::SireVol::TriclinicBox::cubic );
            
            TriclinicBox_exposer.def( 
                "cubic"
                , cubic_function_value
                , ( bp::arg("d") )
                , bp::release_gil_policy()
                , "Return a cubic TriclinicBox with image distance d." );
        
        }
        { //::SireVol::TriclinicBox::gamma
        
            typedef double ( ::SireVol::TriclinicBox::*gamma_function_type)(  ) const;
            gamma_function_type gamma_function_value( &::SireVol::TriclinicBox::gamma );
            
            TriclinicBox_exposer.def( 
                "gamma"
                , gamma_function_value
                , bp::release_gil_policy()
                , "Return the angle between v1 and v0 in degrees." );
        
        }
        { //::SireVol::TriclinicBox::getBoxCenter
        
            typedef ::SireMaths::Vector ( ::SireVol::TriclinicBox::*getBoxCenter_function_type)( ::SireMaths::Vector const & ) const;
            getBoxCenter_function_type getBoxCenter_function_value( &::SireVol::TriclinicBox::getBoxCenter );
            
            TriclinicBox_exposer.def( 
                "getBoxCenter"
                , getBoxCenter_function_value
                , ( bp::arg("p") )
                , bp::release_gil_policy()
                , "Return the center of the box that contains the point p assuming\nthat the center for the central box is located at the origin" );
        
        }
        { //::SireVol::TriclinicBox::getBoxCenter
        
            typedef ::SireMaths::Vector ( ::SireVol::TriclinicBox::*getBoxCenter_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const & ) const;
            getBoxCenter_function_type getBoxCenter_function_value( &::SireVol::TriclinicBox::getBoxCenter );
            
            TriclinicBox_exposer.def( 
                "getBoxCenter"
                , getBoxCenter_function_value
                , ( bp::arg("p"), bp::arg("center") )
                , bp::release_gil_policy()
                , "Return the center of the box that contains the point p assuming\nthat the center for the central box is located at center" );
        
        }
        { //::SireVol::TriclinicBox::getCopiesWithin
        
            typedef ::QList< boost::tuples::tuple< double, SireVol::CoordGroup, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > > ( ::SireVol::TriclinicBox::*getCopiesWithin_function_type)( ::SireVol::CoordGroup const &,::SireVol::CoordGroup const &,double ) const;
            getCopiesWithin_function_type getCopiesWithin_function_value( &::SireVol::TriclinicBox::getCopiesWithin );
            
            TriclinicBox_exposer.def( 
                "getCopiesWithin"
                , getCopiesWithin_function_value
                , ( bp::arg("group"), bp::arg("center"), bp::arg("dist") )
                , bp::release_gil_policy()
                , "Return a list of copies of CoordGroup group that are within\ndistance of the CoordGroup center, translating group so that\nit has the right coordinates to be around center. Note that multiple\ncopies of group may be returned in this is a periodic space and\nthere are multiple periodic replicas of group within dist of\ncenter. The copies of group are returned together with the\nminimum distance between that periodic replica and center.\nIf there are no periodic replicas of group that are within\ndist of center, then an empty list is returned." );
        
        }
        { //::SireVol::TriclinicBox::getImagesWithin
        
            typedef ::QVector< SireMaths::Vector > ( ::SireVol::TriclinicBox::*getImagesWithin_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const &,double ) const;
            getImagesWithin_function_type getImagesWithin_function_value( &::SireVol::TriclinicBox::getImagesWithin );
            
            TriclinicBox_exposer.def( 
                "getImagesWithin"
                , getImagesWithin_function_value
                , ( bp::arg("point"), bp::arg("center"), bp::arg("dist") )
                , bp::release_gil_policy()
                , "Return all periodic images of point with respect to center within\ndist distance of center" );
        
        }
        { //::SireVol::TriclinicBox::getMinimumImage
        
            typedef ::QVector< SireMaths::Vector > ( ::SireVol::TriclinicBox::*getMinimumImage_function_type)( ::QVector< SireMaths::Vector > const &,::SireMaths::Vector const & ) const;
            getMinimumImage_function_type getMinimumImage_function_value( &::SireVol::TriclinicBox::getMinimumImage );
            
            TriclinicBox_exposer.def( 
                "getMinimumImage"
                , getMinimumImage_function_value
                , ( bp::arg("coords"), bp::arg("center") )
                , bp::release_gil_policy()
                , "Return the copy of the point point which is the closest minimum image\nto center" );
        
        }
        { //::SireVol::TriclinicBox::getMinimumImage
        
            typedef ::SireMaths::Vector ( ::SireVol::TriclinicBox::*getMinimumImage_function_type)( ::SireMaths::Vector const &,::SireMaths::Vector const & ) const;
            getMinimumImage_function_type getMinimumImage_function_value( &::SireVol::TriclinicBox::getMinimumImage );
            
            TriclinicBox_exposer.def( 
                "getMinimumImage"
                , getMinimumImage_function_value
                , ( bp::arg("point"), bp::arg("center") )
                , bp::release_gil_policy()
                , "Return the copy of the point point which is the closest minimum image\nto center" );
        
        }
        { //::SireVol::TriclinicBox::getMinimumImage
        
            typedef ::SireVol::CoordGroup ( ::SireVol::TriclinicBox::*getMinimumImage_function_type)( ::SireVol::CoordGroup const &,::SireMaths::Vector const & ) const;
            getMinimumImage_function_type getMinimumImage_function_value( &::SireVol::TriclinicBox::getMinimumImage );
            
            TriclinicBox_exposer.def( 
                "getMinimumImage"
                , getMinimumImage_function_value
                , ( bp::arg("group"), bp::arg("center") )
                , bp::release_gil_policy()
                , "Return the closest periodic copy of group to the point point,\naccording to the minimum image convention. The effect of this is\nto move group into the box which is now centered on point" );
        
        }
        { //::SireVol::TriclinicBox::getMinimumImage
        
            typedef ::SireVol::CoordGroupArray ( ::SireVol::TriclinicBox::*getMinimumImage_function_type)( ::SireVol::CoordGroupArray const &,::SireMaths::Vector const &,bool ) const;
            getMinimumImage_function_type getMinimumImage_function_value( &::SireVol::TriclinicBox::getMinimumImage );
            
            TriclinicBox_exposer.def( 
                "getMinimumImage"
                , getMinimumImage_function_value
                , ( bp::arg("groups"), bp::arg("center"), bp::arg("translate_as_one")=(bool)(false) )
                , "Return the closest periodic copy of each group in groups to the\npoint point, according to the minimum image convention.\nThe effect of this is to move each group into the box which is\nnow centered on point. If translate_as_one is true,\nthen this treats all groups as being part of one larger\ngroup, and so it translates it together. This is useful\nto get the minimum image of a molecule as a whole, rather\nthan breaking the molecule across a box boundary" );
        
        }
        { //::SireVol::TriclinicBox::getMinimumImage
        
            typedef ::SireVol::AABox ( ::SireVol::TriclinicBox::*getMinimumImage_function_type)( ::SireVol::AABox const &,::SireMaths::Vector const & ) const;
            getMinimumImage_function_type getMinimumImage_function_value( &::SireVol::TriclinicBox::getMinimumImage );
            
            TriclinicBox_exposer.def( 
                "getMinimumImage"
                , getMinimumImage_function_value
                , ( bp::arg("aabox"), bp::arg("center") )
                , bp::release_gil_policy()
                , "Return the copy of the triclinic box which is the closest minimum image\nto center" );
        
        }
        { //::SireVol::TriclinicBox::getRandomPoint
        
            typedef ::SireMaths::Vector ( ::SireVol::TriclinicBox::*getRandomPoint_function_type)( ::SireMaths::Vector const &,::SireMaths::RanGenerator const & ) const;
            getRandomPoint_function_type getRandomPoint_function_value( &::SireVol::TriclinicBox::getRandomPoint );
            
            TriclinicBox_exposer.def( 
                "getRandomPoint"
                , getRandomPoint_function_value
                , ( bp::arg("center"), bp::arg("generator") )
                , bp::release_gil_policy()
                , "Return a random point within the box (placing the center of the box\nis at the center center)" );
        
        }
        { //::SireVol::TriclinicBox::isCartesian
        
            typedef bool ( ::SireVol::TriclinicBox::*isCartesian_function_type)(  ) const;
            isCartesian_function_type isCartesian_function_value( &::SireVol::TriclinicBox::isCartesian );
            
            TriclinicBox_exposer.def( 
                "isCartesian"
                , isCartesian_function_value
                , bp::release_gil_policy()
                , "In general, a triclinic box isnt Cartesian." );
        
        }
        { //::SireVol::TriclinicBox::isPeriodic
        
            typedef bool ( ::SireVol::TriclinicBox::*isPeriodic_function_type)(  ) const;
            isPeriodic_function_type isPeriodic_function_value( &::SireVol::TriclinicBox::isPeriodic );
            
            TriclinicBox_exposer.def( 
                "isPeriodic"
                , isPeriodic_function_value
                , bp::release_gil_policy()
                , "A Triclinic box is periodic" );
        
        }
        { //::SireVol::TriclinicBox::isReduced
        
            typedef bool ( ::SireVol::TriclinicBox::*isReduced_function_type)(  ) const;
            isReduced_function_type isReduced_function_value( &::SireVol::TriclinicBox::isReduced );
            
            TriclinicBox_exposer.def( 
                "isReduced"
                , isReduced_function_value
                , bp::release_gil_policy()
                , "Whether an automatic lattice reduction has been performed." );
        
        }
        { //::SireVol::TriclinicBox::isRotated
        
            typedef bool ( ::SireVol::TriclinicBox::*isRotated_function_type)(  ) const;
            isRotated_function_type isRotated_function_value( &::SireVol::TriclinicBox::isRotated );
            
            TriclinicBox_exposer.def( 
                "isRotated"
                , isRotated_function_value
                , bp::release_gil_policy()
                , "Whether the triclinic cell has been rotated to comply with the constraints\nof molecular dynamics engines, i.e. vector0 aligned with x axis, vector1\nin x-y plane, and vector2 with positive z component.\n" );
        
        }
        { //::SireVol::TriclinicBox::minimumDistance
        
            typedef double ( ::SireVol::TriclinicBox::*minimumDistance_function_type)( ::SireVol::CoordGroup const &,::SireVol::CoordGroup const & ) const;
            minimumDistance_function_type minimumDistance_function_value( &::SireVol::TriclinicBox::minimumDistance );
            
            TriclinicBox_exposer.def( 
                "minimumDistance"
                , minimumDistance_function_value
                , ( bp::arg("group0"), bp::arg("group1") )
                , bp::release_gil_policy()
                , "Return the minimum distance between the points in group0 and group1.\nIf this is a periodic space then this uses the minimum image convention\n(i.e. the minimum distance between the closest periodic replicas are\nused)" );
        
        }
        { //::SireVol::TriclinicBox::minimumDistance
        
            typedef double ( ::SireVol::TriclinicBox::*minimumDistance_function_type)( ::SireVol::AABox const &,::SireVol::AABox const & ) const;
            minimumDistance_function_type minimumDistance_function_value( &::SireVol::TriclinicBox::minimumDistance );
            
            TriclinicBox_exposer.def( 
                "minimumDistance"
                , minimumDistance_function_value
                , ( bp::arg("box0"), bp::arg("box1") )
                , bp::release_gil_policy()
                , "Return the minimum distance between the two boxes" );
        
        }
        TriclinicBox_exposer.def( bp::self != bp::self );
        { //::SireVol::TriclinicBox::operator=
        
            typedef ::SireVol::TriclinicBox & ( ::SireVol::TriclinicBox::*assign_function_type)( ::SireVol::TriclinicBox const & ) ;
            assign_function_type assign_function_value( &::SireVol::TriclinicBox::operator= );
            
            TriclinicBox_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        TriclinicBox_exposer.def( bp::self == bp::self );
        { //::SireVol::TriclinicBox::reduce
        
            typedef void ( ::SireVol::TriclinicBox::*reduce_function_type)( double ) ;
            reduce_function_type reduce_function_value( &::SireVol::TriclinicBox::reduce );
            
            TriclinicBox_exposer.def( 
                "reduce"
                , reduce_function_value
                , ( bp::arg("bias")=0. )
                , "Perform a lattice reduction on the triclinic cell.\n\nPar:am bias\nThe bias to use when rounding during the lattice reduction.\nNegative values biases towards left-tilting boxes, whereas\npositive values biases towards right-tilting boxes. This can\nbe used to ensure that rounding is performed in a consistent\ndirection, avoiding oscillation when the TriclinicBox is\ninstantiated from box vectors, or dimensions and angles, that\nhave been read from fixed-precision input files.\n" );
        
        }
        { //::SireVol::TriclinicBox::rhombicDodecahedronHexagon
        
            typedef ::SireVol::TriclinicBox ( *rhombicDodecahedronHexagon_function_type )( double,bool,bool );
            rhombicDodecahedronHexagon_function_type rhombicDodecahedronHexagon_function_value( &::SireVol::TriclinicBox::rhombicDodecahedronHexagon );
            
            TriclinicBox_exposer.def( 
                "rhombicDodecahedronHexagon"
                , rhombicDodecahedronHexagon_function_value
                , ( bp::arg("d"), bp::arg("auto_rotate")=(bool)(true), bp::arg("auto_reduce")=(bool)(true) )
                , "Return a hexagonal rhombic dodecahedron TriclinicBox with image distance d." );
        
        }
        { //::SireVol::TriclinicBox::rhombicDodecahedronSquare
        
            typedef ::SireVol::TriclinicBox ( *rhombicDodecahedronSquare_function_type )( double,bool,bool );
            rhombicDodecahedronSquare_function_type rhombicDodecahedronSquare_function_value( &::SireVol::TriclinicBox::rhombicDodecahedronSquare );
            
            TriclinicBox_exposer.def( 
                "rhombicDodecahedronSquare"
                , rhombicDodecahedronSquare_function_value
                , ( bp::arg("d"), bp::arg("auto_rotate")=(bool)(true), bp::arg("auto_reduce")=(bool)(true) )
                , "Return a square rhombic dodecahedron TriclinicBox with image distance d." );
        
        }
        { //::SireVol::TriclinicBox::rotate
        
            typedef void ( ::SireVol::TriclinicBox::*rotate_function_type)( double ) ;
            rotate_function_type rotate_function_value( &::SireVol::TriclinicBox::rotate );
            
            TriclinicBox_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("precision")=0. )
                , "Rotate the triclinic cell to comply with the constraints of certain\nmolecular dynamics engines, i.e. such that vector0 is aligned with\nthe x axis, vector1, lies in the x-y plane, and vector2 has a positive\nz component.\n\nPar:am precision\nThe precision to use when sorting the lattice vectors based on\ntheir magnitude. This can be used to prevent unwanted rotation\nwhen using input fixed-precision ascii molecular input files.\n" );
        
        }
        { //::SireVol::TriclinicBox::rotationMatrix
        
            typedef ::SireMaths::Matrix const & ( ::SireVol::TriclinicBox::*rotationMatrix_function_type)(  ) const;
            rotationMatrix_function_type rotationMatrix_function_value( &::SireVol::TriclinicBox::rotationMatrix );
            
            TriclinicBox_exposer.def( 
                "rotationMatrix"
                , rotationMatrix_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the rotation matrix." );
        
        }
        { //::SireVol::TriclinicBox::setVolume
        
            typedef ::SireVol::SpacePtr ( ::SireVol::TriclinicBox::*setVolume_function_type)( ::SireUnits::Dimension::Volume ) const;
            setVolume_function_type setVolume_function_value( &::SireVol::TriclinicBox::setVolume );
            
            TriclinicBox_exposer.def( 
                "setVolume"
                , setVolume_function_value
                , ( bp::arg("volume") )
                , bp::release_gil_policy()
                , "Set the volume of the triclinic box." );
        
        }
        { //::SireVol::TriclinicBox::toString
        
            typedef ::QString ( ::SireVol::TriclinicBox::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireVol::TriclinicBox::toString );
            
            TriclinicBox_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this space" );
        
        }
        { //::SireVol::TriclinicBox::truncatedOctahedron
        
            typedef ::SireVol::TriclinicBox ( *truncatedOctahedron_function_type )( double,bool,bool );
            truncatedOctahedron_function_type truncatedOctahedron_function_value( &::SireVol::TriclinicBox::truncatedOctahedron );
            
            TriclinicBox_exposer.def( 
                "truncatedOctahedron"
                , truncatedOctahedron_function_value
                , ( bp::arg("d"), bp::arg("auto_rotate")=(bool)(true), bp::arg("auto_reduce")=(bool)(true) )
                , "Return a truncated octahedron with image distance d." );
        
        }
        { //::SireVol::TriclinicBox::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::TriclinicBox::typeName );
            
            TriclinicBox_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::TriclinicBox::vector0
        
            typedef ::SireMaths::Vector const & ( ::SireVol::TriclinicBox::*vector0_function_type)(  ) const;
            vector0_function_type vector0_function_value( &::SireVol::TriclinicBox::vector0 );
            
            TriclinicBox_exposer.def( 
                "vector0"
                , vector0_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the first box vector." );
        
        }
        { //::SireVol::TriclinicBox::vector1
        
            typedef ::SireMaths::Vector const & ( ::SireVol::TriclinicBox::*vector1_function_type)(  ) const;
            vector1_function_type vector1_function_value( &::SireVol::TriclinicBox::vector1 );
            
            TriclinicBox_exposer.def( 
                "vector1"
                , vector1_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the second box vector." );
        
        }
        { //::SireVol::TriclinicBox::vector2
        
            typedef ::SireMaths::Vector const & ( ::SireVol::TriclinicBox::*vector2_function_type)(  ) const;
            vector2_function_type vector2_function_value( &::SireVol::TriclinicBox::vector2 );
            
            TriclinicBox_exposer.def( 
                "vector2"
                , vector2_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the third box vector." );
        
        }
        { //::SireVol::TriclinicBox::volume
        
            typedef ::SireUnits::Dimension::Volume ( ::SireVol::TriclinicBox::*volume_function_type)(  ) const;
            volume_function_type volume_function_value( &::SireVol::TriclinicBox::volume );
            
            TriclinicBox_exposer.def( 
                "volume"
                , volume_function_value
                , bp::release_gil_policy()
                , "Get the volume of the triclinic box." );
        
        }
        TriclinicBox_exposer.staticmethod( "cubic" );
        TriclinicBox_exposer.staticmethod( "rhombicDodecahedronHexagon" );
        TriclinicBox_exposer.staticmethod( "rhombicDodecahedronSquare" );
        TriclinicBox_exposer.staticmethod( "truncatedOctahedron" );
        TriclinicBox_exposer.staticmethod( "typeName" );
        TriclinicBox_exposer.def( "__copy__", &__copy__);
        TriclinicBox_exposer.def( "__deepcopy__", &__copy__);
        TriclinicBox_exposer.def( "clone", &__copy__);
        TriclinicBox_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireVol::TriclinicBox >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TriclinicBox_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireVol::TriclinicBox >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TriclinicBox_exposer.def_pickle(sire_pickle_suite< ::SireVol::TriclinicBox >());
        TriclinicBox_exposer.def( "__str__", &__str__< ::SireVol::TriclinicBox > );
        TriclinicBox_exposer.def( "__repr__", &__str__< ::SireVol::TriclinicBox > );
    }

}
