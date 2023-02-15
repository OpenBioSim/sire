// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "UnitTest.pypp.hpp"

namespace bp = boost::python;

#include "SireError/exception.h"

#include "unittest.h"

#include <QElapsedTimer>

#include <QTextStream>

#include "unittest.h"

SireBase::UnitTest __copy__(const SireBase::UnitTest &other){ return SireBase::UnitTest(other); }

const char* pvt_get_name(const SireBase::UnitTest&){ return "SireBase::UnitTest";}

#include "Helpers/release_gil_policy.hpp"

void register_UnitTest_class(){

    { //::SireBase::UnitTest
        typedef bp::class_< SireBase::UnitTest > UnitTest_exposer_t;
        UnitTest_exposer_t UnitTest_exposer = UnitTest_exposer_t( "UnitTest", "This class is used to register a unit test\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope UnitTest_scope( UnitTest_exposer );
        UnitTest_exposer.def( bp::init< SireBase::UnitTest const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireBase::UnitTest::errorString

            typedef ::QString ( ::SireBase::UnitTest::*errorString_function_type)(  ) ;
            errorString_function_type errorString_function_value( &::SireBase::UnitTest::errorString );

            UnitTest_exposer.def(
                "errorString"
                , errorString_function_value
                , bp::release_gil_policy()
                , "Return the error string" );

        }
        { //::SireBase::UnitTest::name

            typedef ::QString ( ::SireBase::UnitTest::*name_function_type)(  ) const;
            name_function_type name_function_value( &::SireBase::UnitTest::name );

            UnitTest_exposer.def(
                "name"
                , name_function_value
                , bp::release_gil_policy()
                , "Return the name of the test" );

        }
        { //::SireBase::UnitTest::run

            typedef bool ( ::SireBase::UnitTest::*run_function_type)( int,bool ) ;
            run_function_type run_function_value( &::SireBase::UnitTest::run );

            UnitTest_exposer.def(
                "run"
                , run_function_value
                , ( bp::arg("nrepeats")=(int)(1), bp::arg("verbose")=(bool)(false) )
                , "Run this test a total of nrepeats times, printing out information\nif verbose is true" );

        }
        { //::SireBase::UnitTest::runAll

            typedef int ( *runAll_function_type )( bool );
            runAll_function_type runAll_function_value( &::SireBase::UnitTest::runAll );

            UnitTest_exposer.def(
                "runAll"
                , runAll_function_value
                , ( bp::arg("verbose")=(bool)(false) )
                , "Run all of the tests one after another, printing out the results to the\nscreen and returning the number of tests that failed" );

        }
        { //::SireBase::UnitTest::runTime

            typedef ::quint64 ( ::SireBase::UnitTest::*runTime_function_type)(  ) ;
            runTime_function_type runTime_function_value( &::SireBase::UnitTest::runTime );

            UnitTest_exposer.def(
                "runTime"
                , runTime_function_value
                , bp::release_gil_policy()
                , "Return the time taken to run the test, in nanoseconds" );

        }
        { //::SireBase::UnitTest::tests

            typedef ::QList< std::shared_ptr< SireBase::UnitTest > > ( *tests_function_type )(  );
            tests_function_type tests_function_value( &::SireBase::UnitTest::tests );

            UnitTest_exposer.def(
                "tests"
                , tests_function_value
                , bp::release_gil_policy()
                , "Return all of the tests that have been registered" );

        }
        { //::SireBase::UnitTest::wasError

            typedef bool ( ::SireBase::UnitTest::*wasError_function_type)(  ) ;
            wasError_function_type wasError_function_value( &::SireBase::UnitTest::wasError );

            UnitTest_exposer.def(
                "wasError"
                , wasError_function_value
                , bp::release_gil_policy()
                , "Return whether or not the test was an error" );

        }
        { //::SireBase::UnitTest::wasSuccessful

            typedef bool ( ::SireBase::UnitTest::*wasSuccessful_function_type)(  ) ;
            wasSuccessful_function_type wasSuccessful_function_value( &::SireBase::UnitTest::wasSuccessful );

            UnitTest_exposer.def(
                "wasSuccessful"
                , wasSuccessful_function_value
                , bp::release_gil_policy()
                , "Return whether or not the test was successful" );

        }
        UnitTest_exposer.staticmethod( "runAll" );
        UnitTest_exposer.staticmethod( "tests" );
        bp::register_ptr_to_python< std::shared_ptr< SireBase::UnitTest > >();
        UnitTest_exposer.def( "__copy__", &__copy__);
        UnitTest_exposer.def( "__deepcopy__", &__copy__);
        UnitTest_exposer.def( "clone", &__copy__);
        UnitTest_exposer.def( "__str__", &pvt_get_name);
        UnitTest_exposer.def( "__repr__", &pvt_get_name);
    }

}
