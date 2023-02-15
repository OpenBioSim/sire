#ifndef SIREBASE_CONVERT_PROPERTY_HPP
#define SIREBASE_CONVERT_PROPERTY_HPP

// boost/type_traits depends on a depracated function
// Disabling the warning so it doesn't pollute our build log
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-builtins"
#include <boost/type_traits.hpp>
#pragma clang diagnostic pop

#include "generalunitproperty.h"
#include "property.h"
#include "variantproperty.h"

namespace SireBase
{

    namespace detail
    {
        template <int T>
        struct convert_property
        {
            template <class V>
            static SireBase::PropertyPtr convert(const V &value);
        };

        template <>
        struct convert_property<true>
        {
            template <class V>
            static SireBase::PropertyPtr convert(const V &value)
            {
                return SireBase::PropertyPtr(value);
            }
        };

        template <>
        struct convert_property<false>
        {
            template <class V>
            static SireBase::PropertyPtr convert(const V &value)
            {
                return SireBase::PropertyPtr(SireBase::VariantProperty(QVariant::fromValue(value)));
            }
        };

        template <int T>
        struct convert_unit_property
        {
            template <class V>
            static SireBase::PropertyPtr convert(const V &value);
        };

        template <>
        struct convert_unit_property<true>
        {
            template <class V>
            static SireBase::PropertyPtr convert(const V &value)
            {
                return SireBase::PropertyPtr(GeneralUnitProperty(value));
            }
        };

        template <>
        struct convert_unit_property<false>
        {
            template <class V>
            static SireBase::PropertyPtr convert(const V &value)
            {
                return detail::convert_property<boost::is_base_of<SireBase::Property, V>::value>::convert(value);
            }
        };
    } // namespace detail

    template <class T>
    inline SireBase::PropertyPtr convert_property(const T &value)
    {
        return detail::convert_unit_property<boost::is_base_of<SireUnits::Dimension::Unit, T>::value>::convert(value);
    }

    template <>
    inline SireBase::PropertyPtr convert_property(const SireUnits::Dimension::GeneralUnit &value)
    {
        return SireBase::PropertyPtr(GeneralUnitProperty(value));
    }

} // namespace SireBase

#endif
