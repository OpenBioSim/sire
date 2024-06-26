// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Process.pypp.hpp"

namespace bp = boost::python;

#include "sire_process.h"

#include "sire_process.h"

SireBase::Process __copy__(const SireBase::Process &other){ return SireBase::Process(other); }

#include "Helpers/copy.hpp"

const char* pvt_get_name(const SireBase::Process&){ return "SireBase::Process";}

#include "Helpers/release_gil_policy.hpp"

void register_Process_class(){

    { //::SireBase::Process
        typedef bp::class_< SireBase::Process > Process_exposer_t;
        Process_exposer_t Process_exposer = Process_exposer_t( "Process", "This class provides a means to run an external process\n(executable). This provides the equivalent of\nstd::system, but with added error handling and\nsignal handling (which ensures that any child processes\nare killed when Sire exits)\n\nAuthor: Christopher Woods, Lester Hedges\n", bp::init< >("Null constructor") );
        bp::scope Process_scope( Process_exposer );
        Process_exposer.def( bp::init< SireBase::Process const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireBase::Process::hasFinished
        
            typedef bool ( ::SireBase::Process::*hasFinished_function_type)(  ) ;
            hasFinished_function_type hasFinished_function_value( &::SireBase::Process::hasFinished );
            
            Process_exposer.def( 
                "hasFinished"
                , hasFinished_function_value
                , bp::release_gil_policy()
                , "Return whether or not this process has finished running" );
        
        }
        { //::SireBase::Process::isError
        
            typedef bool ( ::SireBase::Process::*isError_function_type)(  ) ;
            isError_function_type isError_function_value( &::SireBase::Process::isError );
            
            Process_exposer.def( 
                "isError"
                , isError_function_value
                , bp::release_gil_policy()
                , "Return whether or not the process exited in error" );
        
        }
        { //::SireBase::Process::isRunning
        
            typedef bool ( ::SireBase::Process::*isRunning_function_type)(  ) ;
            isRunning_function_type isRunning_function_value( &::SireBase::Process::isRunning );
            
            Process_exposer.def( 
                "isRunning"
                , isRunning_function_value
                , bp::release_gil_policy()
                , "Return whether or not the job is running" );
        
        }
        { //::SireBase::Process::kill
        
            typedef void ( ::SireBase::Process::*kill_function_type)(  ) ;
            kill_function_type kill_function_value( &::SireBase::Process::kill );
            
            Process_exposer.def( 
                "kill"
                , kill_function_value
                , bp::release_gil_policy()
                , "Kill this process" );
        
        }
        { //::SireBase::Process::killAll
        
            typedef void ( *killAll_function_type )(  );
            killAll_function_type killAll_function_value( &::SireBase::Process::killAll );
            
            Process_exposer.def( 
                "killAll"
                , killAll_function_value
                , bp::release_gil_policy()
                , "Use this function to kill all of the jobs that are currently running" );
        
        }
        Process_exposer.def( bp::self != bp::self );
        { //::SireBase::Process::operator=
        
            typedef ::SireBase::Process & ( ::SireBase::Process::*assign_function_type)( ::SireBase::Process const & ) ;
            assign_function_type assign_function_value( &::SireBase::Process::operator= );
            
            Process_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Process_exposer.def( bp::self == bp::self );
        { //::SireBase::Process::run
        
            typedef ::SireBase::Process ( *run_function_type )( ::QString const & );
            run_function_type run_function_value( &::SireBase::Process::run );
            
            Process_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("command") )
                , bp::release_gil_policy()
                , "Run the command command and return a Process object that can be\nused to monitor the command" );
        
        }
        { //::SireBase::Process::run
        
            typedef ::SireBase::Process ( *run_function_type )( ::QString const &,::QString const &,::QString const & );
            run_function_type run_function_value( &::SireBase::Process::run );
            
            Process_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("command"), bp::arg("stdout_file"), bp::arg("stderr_file") )
                , bp::release_gil_policy()
                , "Run the command command and return a Process object that can be\nused to monitor the command. Stdout and stderr of the running\nprocess are redirected to the user specified files." );
        
        }
        { //::SireBase::Process::run
        
            typedef ::SireBase::Process ( *run_function_type )( ::QString const &,::QString const & );
            run_function_type run_function_value( &::SireBase::Process::run );
            
            Process_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("command"), bp::arg("arg") )
                , bp::release_gil_policy()
                , "Run the command command with the solitary argument arg" );
        
        }
        { //::SireBase::Process::run
        
            typedef ::SireBase::Process ( *run_function_type )( ::QString const &,::QString const &,::QString const &,::QString const & );
            run_function_type run_function_value( &::SireBase::Process::run );
            
            Process_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("command"), bp::arg("arg"), bp::arg("stdout_file"), bp::arg("stderr_file") )
                , bp::release_gil_policy()
                , "Run the command command with the solitary argument arg.\nStdout and stderr of the running process are redirected to\nthe user specified files." );
        
        }
        { //::SireBase::Process::run
        
            typedef ::SireBase::Process ( *run_function_type )( ::QString const &,::QStringList const & );
            run_function_type run_function_value( &::SireBase::Process::run );
            
            Process_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("command"), bp::arg("arguments") )
                , bp::release_gil_policy()
                , "Run the command command with the arguments arguments, and\nreturn a Process object that can be used to query and control the\njob" );
        
        }
        { //::SireBase::Process::run
        
            typedef ::SireBase::Process ( *run_function_type )( ::QString const &,::QStringList const &,::QString const &,::QString const & );
            run_function_type run_function_value( &::SireBase::Process::run );
            
            Process_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("command"), bp::arg("arguments"), bp::arg("stdout_file"), bp::arg("stderr_file") )
                , bp::release_gil_policy()
                , "Run the command command with the arguments arguments, and\nreturn a Process object that can be used to query and control the\njob. Stdout and stderr of the running process are redirected to\nthe user specified files." );
        
        }
        { //::SireBase::Process::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::Process::typeName );
            
            Process_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Process::wait
        
            typedef void ( ::SireBase::Process::*wait_function_type)(  ) ;
            wait_function_type wait_function_value( &::SireBase::Process::wait );
            
            Process_exposer.def( 
                "wait"
                , wait_function_value
                , bp::release_gil_policy()
                , "Wait until the process has finished" );
        
        }
        { //::SireBase::Process::wait
        
            typedef bool ( ::SireBase::Process::*wait_function_type)( int ) ;
            wait_function_type wait_function_value( &::SireBase::Process::wait );
            
            Process_exposer.def( 
                "wait"
                , wait_function_value
                , ( bp::arg("ms") )
                , bp::release_gil_policy()
                , "Wait until the process has finished, or until ms milliseconds have passed.\nThis returns whether or not the process has finished" );
        
        }
        { //::SireBase::Process::wasKilled
        
            typedef bool ( ::SireBase::Process::*wasKilled_function_type)(  ) ;
            wasKilled_function_type wasKilled_function_value( &::SireBase::Process::wasKilled );
            
            Process_exposer.def( 
                "wasKilled"
                , wasKilled_function_value
                , bp::release_gil_policy()
                , "Return whether or not the process was killed" );
        
        }
        { //::SireBase::Process::what
        
            typedef char const * ( ::SireBase::Process::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireBase::Process::what );
            
            Process_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Process_exposer.staticmethod( "killAll" );
        Process_exposer.staticmethod( "run" );
        Process_exposer.staticmethod( "typeName" );
        Process_exposer.def( "__copy__", &__copy__<SireBase::Process>);
        Process_exposer.def( "__deepcopy__", &__copy__<SireBase::Process>);
        Process_exposer.def( "clone", &__copy__<SireBase::Process>);
        Process_exposer.def( "__str__", &pvt_get_name);
        Process_exposer.def( "__repr__", &pvt_get_name);
    }

}
