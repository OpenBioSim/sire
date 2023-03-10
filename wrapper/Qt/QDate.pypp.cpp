// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include <QString>
#include <QByteArray>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QDateTime>
#include <QLocale>
#include <QUuid>
#include <qnamespace.h>
#include <QVariant>
#include <QUrl>
#include <QBitArray>
#include "QDate.pypp.hpp"

namespace bp = boost::python;

QDate __copy__(const QDate &other){ return QDate(other); }

#include "Qt/qdatastream.hpp"

void register_QDate_class(){

    { //::QDate
        typedef bp::class_< QDate > QDate_exposer_t;
        QDate_exposer_t QDate_exposer = QDate_exposer_t( "QDate", bp::init< >() );
        bp::scope QDate_scope( QDate_exposer );
        QDate_exposer.def( bp::init< int, int, int >(( bp::arg("y"), bp::arg("m"), bp::arg("d") )) );
        { //::QDate::addDays

            typedef ::QDate ( ::QDate::*addDays_function_type )( qint64 ) const;
            addDays_function_type addDays_function_value( &::QDate::addDays );

            QDate_exposer.def(
                "addDays"
                , addDays_function_value
                , ( bp::arg("days") ) );

        }
        { //::QDate::addMonths

            typedef ::QDate ( ::QDate::*addMonths_function_type )( int ) const;
            addMonths_function_type addMonths_function_value( &::QDate::addMonths );

            QDate_exposer.def(
                "addMonths"
                , addMonths_function_value
                , ( bp::arg("months") ) );

        }
        { //::QDate::addYears

            typedef ::QDate ( ::QDate::*addYears_function_type )( int ) const;
            addYears_function_type addYears_function_value( &::QDate::addYears );

            QDate_exposer.def(
                "addYears"
                , addYears_function_value
                , ( bp::arg("years") ) );

        }
        { //::QDate::currentDate

            typedef ::QDate ( *currentDate_function_type )(  );
            currentDate_function_type currentDate_function_value( &::QDate::currentDate );

            QDate_exposer.def(
                "currentDate"
                , currentDate_function_value );

        }
        { //::QDate::day

            typedef int ( ::QDate::*day_function_type )(  ) const;
            day_function_type day_function_value( &::QDate::day );

            QDate_exposer.def(
                "day"
                , day_function_value );

        }
        { //::QDate::dayOfWeek

            typedef int ( ::QDate::*dayOfWeek_function_type )(  ) const;
            dayOfWeek_function_type dayOfWeek_function_value( &::QDate::dayOfWeek );

            QDate_exposer.def(
                "dayOfWeek"
                , dayOfWeek_function_value );

        }
        { //::QDate::dayOfYear

            typedef int ( ::QDate::*dayOfYear_function_type )(  ) const;
            dayOfYear_function_type dayOfYear_function_value( &::QDate::dayOfYear );

            QDate_exposer.def(
                "dayOfYear"
                , dayOfYear_function_value );

        }
        { //::QDate::daysInMonth

            typedef int ( ::QDate::*daysInMonth_function_type )(  ) const;
            daysInMonth_function_type daysInMonth_function_value( &::QDate::daysInMonth );

            QDate_exposer.def(
                "daysInMonth"
                , daysInMonth_function_value );

        }
        { //::QDate::daysInYear

            typedef int ( ::QDate::*daysInYear_function_type )(  ) const;
            daysInYear_function_type daysInYear_function_value( &::QDate::daysInYear );

            QDate_exposer.def(
                "daysInYear"
                , daysInYear_function_value );

        }
        { //::QDate::daysTo

            typedef qint64 ( ::QDate::*daysTo_function_type )( ::QDate const & ) const;
            daysTo_function_type daysTo_function_value( &::QDate::daysTo );

            QDate_exposer.def(
                "daysTo"
                , daysTo_function_value
                , ( bp::arg("arg0") ) );

        }
        { //::QDate::fromJulianDay

            typedef ::QDate ( *fromJulianDay_function_type )( qint64 );
            fromJulianDay_function_type fromJulianDay_function_value( &::QDate::fromJulianDay );

            QDate_exposer.def(
                "fromJulianDay"
                , fromJulianDay_function_value
                , ( bp::arg("jd") ) );

        }
        { //::QDate::fromString

            typedef ::QDate ( *fromString_function_type )( ::QString const &,::QString const & );
            fromString_function_type fromString_function_value( &::QDate::fromString );

            QDate_exposer.def(
                "fromString"
                , fromString_function_value
                , ( bp::arg("s"), bp::arg("format") ) );

        }
        { //::QDate::isLeapYear

            typedef bool ( *isLeapYear_function_type )( int );
            isLeapYear_function_type isLeapYear_function_value( &::QDate::isLeapYear );

            QDate_exposer.def(
                "isLeapYear"
                , isLeapYear_function_value
                , ( bp::arg("year") ) );

        }
        { //::QDate::isNull

            typedef bool ( ::QDate::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::QDate::isNull );

            QDate_exposer.def(
                "isNull"
                , isNull_function_value );

        }
        { //::QDate::isValid

            typedef bool ( ::QDate::*isValid_function_type )(  ) const;
            isValid_function_type isValid_function_value( &::QDate::isValid );

            QDate_exposer.def(
                "isValid"
                , isValid_function_value );

        }
        { //::QDate::isValid

            typedef bool ( *isValid_function_type )( int,int,int );
            isValid_function_type isValid_function_value( &::QDate::isValid );

            QDate_exposer.def(
                "isValid"
                , isValid_function_value
                , ( bp::arg("y"), bp::arg("m"), bp::arg("d") ) );

        }
        { //::QDate::longDayName

            typedef ::QString ( *longDayName_function_type )( int, QDate::MonthNameType );
            longDayName_function_type longDayName_function_value( &::QDate::longDayName );

            QDate_exposer.def(
                "longDayName"
                , longDayName_function_value
                , ( bp::arg("weekday") ) );

        }
        { //::QDate::longMonthName

            typedef ::QString ( *longMonthName_function_type )( int, QDate::MonthNameType );
            longMonthName_function_type longMonthName_function_value( &::QDate::longMonthName );

            QDate_exposer.def(
                "longMonthName"
                , longMonthName_function_value
                , ( bp::arg("month") ) );

        }
        { //::QDate::month

            typedef int ( ::QDate::*month_function_type )(  ) const;
            month_function_type month_function_value( &::QDate::month );

            QDate_exposer.def(
                "month"
                , month_function_value );

        }
        QDate_exposer.def( bp::self != bp::self );
        QDate_exposer.def( bp::self < bp::self );
        QDate_exposer.def( bp::self <= bp::self );
        QDate_exposer.def( bp::self == bp::self );
        QDate_exposer.def( bp::self > bp::self );
        QDate_exposer.def( bp::self >= bp::self );
        { //::QDate::setDate

            typedef bool ( ::QDate::*setDate_function_type )( int,int,int ) ;
            setDate_function_type setDate_function_value( &::QDate::setDate );

            QDate_exposer.def(
                "setDate"
                , setDate_function_value
                , ( bp::arg("year"), bp::arg("month"), bp::arg("date") ) );

        }
        { //::QDate::shortDayName

            typedef ::QString ( *shortDayName_function_type )( int, QDate::MonthNameType );
            shortDayName_function_type shortDayName_function_value( &::QDate::shortDayName );

            QDate_exposer.def(
                "shortDayName"
                , shortDayName_function_value
                , ( bp::arg("weekday") ) );

        }
        { //::QDate::shortMonthName

            typedef ::QString ( *shortMonthName_function_type )( int, QDate::MonthNameType );
            shortMonthName_function_type shortMonthName_function_value( &::QDate::shortMonthName );

            QDate_exposer.def(
                "shortMonthName"
                , shortMonthName_function_value
                , ( bp::arg("month") ) );

        }
        { //::QDate::toJulianDay

            typedef qint64 ( ::QDate::*toJulianDay_function_type )(  ) const;
            toJulianDay_function_type toJulianDay_function_value( &::QDate::toJulianDay );

            QDate_exposer.def(
                "toJulianDay"
                , toJulianDay_function_value );

        }
        { //::QDate::toString

            typedef ::QString ( ::QDate::*toString_function_type )( ::QString const & ) const;
            toString_function_type toString_function_value( &::QDate::toString );

            QDate_exposer.def(
                "toString"
                , toString_function_value
                , ( bp::arg("format") ) );

        }
        { //::QDate::weekNumber

            typedef int ( ::QDate::*weekNumber_function_type )( int * ) const;
            weekNumber_function_type weekNumber_function_value( &::QDate::weekNumber );

            QDate_exposer.def(
                "weekNumber"
                , weekNumber_function_value
                , ( bp::arg("yearNum")=bp::object() ) );

        }
        { //::QDate::year

            typedef int ( ::QDate::*year_function_type )(  ) const;
            year_function_type year_function_value( &::QDate::year );

            QDate_exposer.def(
                "year"
                , year_function_value );

        }
        QDate_exposer.staticmethod( "currentDate" );
        QDate_exposer.staticmethod( "fromJulianDay" );
        QDate_exposer.staticmethod( "fromString" );
        QDate_exposer.staticmethod( "isLeapYear" );
        QDate_exposer.staticmethod( "isValid" );
        QDate_exposer.staticmethod( "longDayName" );
        QDate_exposer.staticmethod( "longMonthName" );
        QDate_exposer.staticmethod( "shortDayName" );
        QDate_exposer.staticmethod( "shortMonthName" );
        QDate_exposer.def( "__copy__", &__copy__);
        QDate_exposer.def( "__deepcopy__", &__copy__);
        QDate_exposer.def( "clone", &__copy__);
        QDate_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::QDate >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QDate_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::QDate >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
    }

}
