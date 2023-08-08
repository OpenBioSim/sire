/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include <QHash>

#include "properties.h"
#include "refcountdata.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <memory>

#include <QDebug>

using namespace SireBase;
using namespace SireStream;

namespace SireBase
{

    namespace detail
    {

        /** This is the data for a Properties object */
        class PropertiesData : public RefCountData
        {
        public:
            PropertiesData();
            PropertiesData(const PropertiesData &other);

            ~PropertiesData();

            PropertiesData &operator=(const PropertiesData &other);

            bool operator==(const PropertiesData &other) const;
            bool operator!=(const PropertiesData &other) const;

            bool hasMetadata() const;
            bool hasPropsMetadata() const;
            bool hasLinks() const;

            /** The metadata for this set of properties itself */
            std::shared_ptr<Properties> metadata;

            /** The collection of properties, indexed by name */
            QHash<QString, PropertyPtr> properties;

            /** The metadata for each property, indexed by name */
            std::shared_ptr<QHash<QString, Properties>> props_metadata;

            /** The set of linked properties */
            std::shared_ptr<QHash<QString, QString>> prop_links;
        };

    } // namespace detail

} // namespace SireBase

using namespace SireBase::detail;

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const PropertiesData &props)
{
    SharedDataStream sds(ds);

    if (props.metadata.get() == 0)
    {
        sds << Properties();
    }
    else
    {
        sds << *(props.metadata);
    }

    sds << quint32(props.properties.count());

    for (QHash<QString, PropertyPtr>::const_iterator it = props.properties.constBegin();
         it != props.properties.constEnd(); ++it)
    {
        sds << it.key() << it.value();
    }

    if (props.props_metadata.get() == 0)
    {
        sds << quint32(0);
    }
    else
    {
        sds << quint32(props.props_metadata->count());

        for (QHash<QString, Properties>::const_iterator it = props.props_metadata->constBegin();
             it != props.props_metadata->constEnd(); ++it)
        {
            sds << it.key() << it.value();
        }
    }

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, PropertiesData &props)
{
    SharedDataStream sds(ds);

    Properties metadata;

    sds >> metadata;

    if (metadata.isEmpty())
    {
        props.metadata.reset();
    }
    else
    {
        props.metadata.reset(new Properties(metadata));
    }

    quint32 nprops;

    sds >> nprops;

    props.properties.clear();

    if (nprops > 0)
    {
        props.properties.reserve(nprops);

        for (quint32 i = 0; i < nprops; ++i)
        {
            QString key;
            PropertyPtr value;

            sds >> key >> value;

            props.properties.insert(key, value);
        }
    }

    sds >> nprops;

    if (nprops <= 0)
    {
        props.props_metadata.reset();
    }
    else
    {
        props.props_metadata.reset(new QHash<QString, Properties>());
    }

    if (nprops > 0)
    {
        props.props_metadata->reserve(nprops);

        for (quint32 i = 0; i < nprops; ++i)
        {
            QString key;
            Properties value;

            sds >> key >> value;

            props.props_metadata->insert(key, value);
        }
    }

    return ds;
}

/** Null constructor */
PropertiesData::PropertiesData()
{
}

/** Copy constructor */
PropertiesData::PropertiesData(const PropertiesData &other)
{
    this->operator=(other);
}

/** Destructor */
PropertiesData::~PropertiesData()
{
}

/** Assignment operator */
PropertiesData &PropertiesData::operator=(const PropertiesData &other)
{
    if (this != &other)
    {
        properties = other.properties;

        if (other.metadata.get() != 0)
        {
            metadata.reset(new Properties(*(other.metadata)));
        }

        if (other.props_metadata.get() != 0)
        {
            props_metadata.reset(new QHash<QString, Properties>(*(other.props_metadata)));
        }

        if (other.prop_links.get() != 0)
        {
            prop_links.reset(new QHash<QString, QString>(*(other.prop_links)));
        }
    }

    return *this;
}

bool PropertiesData::hasMetadata() const
{
    if (metadata.get() != 0)
        return not metadata->isEmpty();
    else
        return false;
}

bool PropertiesData::hasPropsMetadata() const
{
    if (props_metadata.get() != 0)
        return not props_metadata->isEmpty();
    else
        return false;
}

bool PropertiesData::hasLinks() const
{
    if (prop_links.get() != 0)
        return not prop_links->isEmpty();
    else
        return false;
}

/** Comparison operator */
bool PropertiesData::operator==(const PropertiesData &other) const
{
    if (this == &other)
        return true;

    if (properties == other.properties)
    {
        if (this->hasMetadata())
        {
            if (other.hasMetadata())
            {
                if (*metadata != *(other.metadata))
                    return false;
            }
            else
                return false;
        }
        else if (other.hasMetadata())
            return false;

        if (this->hasPropsMetadata())
        {
            if (other.hasPropsMetadata())
            {
                if (*props_metadata != *(other.props_metadata))
                    return false;
            }
            else
                return false;
        }
        else if (other.hasPropsMetadata())
            return false;

        if (this->hasLinks())
        {
            if (other.hasLinks())
            {
                if (*prop_links != *(other.prop_links))
                    return false;
            }
            else
                return false;
        }
        else if (other.hasLinks())
            return false;

        return true;
    }
    else
        return false;
}

/** Comparison operator */
bool PropertiesData::operator!=(const PropertiesData &other) const
{
    return not this->operator==(other);
}

/////////////
///////////// Implementation of Properties
/////////////

static const RegisterMetaType<Properties> r_props;

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const Properties &props)
{
    const bool has_links = props.hasLinks();

    if (has_links)
        writeHeader(ds, r_props, 2);
    else
        writeHeader(ds, r_props, 1);

    SharedDataStream sds(ds);

    if (props.isEmpty())
    {
        // just serialise the number 0 to indicate that this is an empty
        // properties object - this is to prevent circular references
        sds << qint32(0);
    }
    else if (has_links)
    {
        sds << qint32(1) << props.d
            << *(props.d->prop_links)
            << static_cast<const Property &>(props);
    }
    else
    {
        sds << qint32(1) << props.d
            << static_cast<const Property &>(props);
    }

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, Properties &props)
{
    VersionID v = readHeader(ds, r_props);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        qint32 not_empty;

        sds >> not_empty;

        if (not_empty)
        {
            QHash<QString, QString> links;

            sds >> props.d >> links >> static_cast<Property &>(props);

            if (not props.isEmpty())
                props.d->prop_links.reset(new QHash<QString, QString>(links));
        }
        else
        {
            // this is an empty Properties object
            props = Properties();
        }
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        qint32 not_empty;

        sds >> not_empty;

        if (not_empty)
        {
            sds >> props.d >> static_cast<Property &>(props);
        }
        else
        {
            // this is an empty Properties object
            props = Properties();
        }
    }
    else
        throw version_error(v, "1, 2", r_props, CODELOC);

    return ds;
}

/** Private constructor - not used any more */
Properties::Properties(bool) : ConcreteProperty<Properties, Property>(), d(0)
{
}

/** Null constructor - construct an empty set of properties */
Properties::Properties() : ConcreteProperty<Properties, Property>(), d(0)
{
}

/** Copy constructor */
Properties::Properties(const Properties &other) : ConcreteProperty<Properties, Property>(other), d(other.d)
{
}

/** Destructor */
Properties::~Properties()
{
}

/** Assignment operator */
Properties &Properties::operator=(const Properties &other)
{
    d = other.d;
    Property::operator=(other);

    return *this;
}

/** Comparison operator */
bool Properties::operator==(const Properties &other) const
{
    if (this->isEmpty() or other.isEmpty())
    {
        return d == other.d;
    }
    else
    {
        return d == other.d or *d == *(other.d);
    }
}

/** Comparison operator */
bool Properties::operator!=(const Properties &other) const
{
    return not this->operator==(other);
}

/** Return whether this is empty (has no values) */
bool Properties::isEmpty() const
{
    return d.constData() == 0;
}

/** Add a link from the property 'key' to the property 'linked_property'.
 *  The linked_property will be returned if there is no property
 *  called 'key' in this set.
 *
 *  Note that the linked property must already be contained in this set.
 */
void Properties::addLink(const QString &key, const QString &linked_property)
{
    this->assertContainsProperty(linked_property);

    if (d->properties.contains(key))
        throw SireError::invalid_key(QObject::tr(
                                         "You cannot create the link '%1' as a property with this name "
                                         "already exists.")
                                         .arg(key),
                                     CODELOC);

    if (d->prop_links.get() == 0)
        d->prop_links.reset(new QHash<QString, QString>());

    d->prop_links->insert(key, linked_property);
}

/** Remove the link associated with the key 'key' */
void Properties::removeLink(const QString &key)
{
    if (not this->isEmpty())
    {
        if (d->prop_links.get() == 0)
            return;

        d->prop_links->remove(key);

        if (d->prop_links->isEmpty())
            d->prop_links.reset();
    }
}

/** Remove all property links from this set */
void Properties::removeAllLinks()
{
    if (not this->isEmpty())
    {
        d->prop_links.reset();
    }
}

/** Return whether or not there are any property links */
bool Properties::hasLinks() const
{
    if (this->isEmpty())
        return false;
    else
        return d->prop_links.get() != 0;
}

/** Return all of the property links */
QHash<QString, QString> Properties::getLinks() const
{
    if (not this->isEmpty())
    {
        if (d->prop_links.get() != 0)
            return *(d->prop_links);
    }

    return QHash<QString, QString>();
}

/** Return the keys for all of the properties in this set */
QStringList Properties::propertyKeys() const
{
    if (this->isEmpty())
        return QStringList();
    else
        return d->properties.keys();
}

/** Return an iterator pointing to the first property in this set */
Properties::const_iterator Properties::begin() const
{
    if (this->isEmpty())
        return Properties::const_iterator();
    else
        return d->properties.begin();
}

/** Return an iterator pointing to the first property in this set */
Properties::const_iterator Properties::constBegin() const
{
    if (this->isEmpty())
        return Properties::const_iterator();
    else
        return d->properties.constBegin();
}

/** Internal function used to get the linked property name
 *  for 'key'. This returns an empty string if there
 *  is no link
 */
QString Properties::_getLink(const QString &key) const
{
    if (this->isEmpty())
        return QString();
    else if (d->prop_links.get() == 0)
        return QString();
    else
        return d->prop_links->value(key, QString());
}

/** Return an iterator pointing to the property with key 'key', or
    Properties::end() if there is no such property */
Properties::const_iterator Properties::find(const QString &key) const
{
    if (this->isEmpty())
        return Properties::const_iterator();
    else
    {
        auto it = d->properties.find(key);

        if (it == d->properties.constEnd())
        {
            auto linked = this->_getLink(key);

            if (not linked.isEmpty())
                it = d->properties.find(linked);
        }

        return it;
    }
}

/** Return an iterator pointing to the property with key 'key', or
    Properties::end() if there is no such property */
Properties::const_iterator Properties::constFind(const QString &key) const
{
    if (this->isEmpty())
        return Properties::const_iterator();
    else
        return d->properties.constFind(key);
}

/** Return an iterator pointing one beyond the last property in this set */
Properties::const_iterator Properties::end() const
{
    if (this->isEmpty())
        return Properties::const_iterator();
    else
        return d->properties.end();
}

/** Return an iterator pointing one beyond the last property in this set */
Properties::const_iterator Properties::constEnd() const
{
    if (this->isEmpty())
        return Properties::const_iterator();
    else
        return d->properties.end();
}

/** Assert that this set contains a property with key 'key'

    \throw SireBase::missing_property
*/
void Properties::assertContainsProperty(const PropertyName &key) const
{
    if (this->isEmpty())
    {
        throw SireBase::missing_property(QObject::tr("There is no property with key \"%1\" as this set of properties is empty.")
                                             .arg(key.source()),
                                         CODELOC);
    }
    else if (key.hasSource() and not d->properties.contains(key.source()) and not key.hasDefaultValue())
    {
        auto linked = this->_getLink(key.source());

        if (linked.isEmpty() or not d->properties.contains(linked))
            throw SireBase::missing_property(QObject::tr("There is no property with key \"%1\". Available keys are ( %2 ).")
                                                 .arg(key.source(), this->propertyKeys().join(", ")),
                                             CODELOC);
    }
}

/** Return the list of metadata keys */
QStringList Properties::metadataKeys() const
{
    if (this->isEmpty())
        return QStringList();
    else if (d->metadata.get() == 0)
        return QStringList();
    else
        return d->metadata->propertyKeys();
}

/** Assert that this contains the metadata at metakey 'metakey' */
void Properties::assertContainsMetadata(const PropertyName &metakey) const
{
    if (this->isEmpty())
    {
        throw SireBase::missing_property(QObject::tr("There is no metadata with metakey \"%1\" as this set of properties is empty.")
                                             .arg(metakey.source()),
                                         CODELOC);
    }
    else if (d->metadata.get() == 0)
    {
        throw SireBase::missing_property(QObject::tr("There is no metadata with metakey \"%1\" as there is no metadata.")
                                             .arg(metakey.source()),
                                         CODELOC);
    }
    else
    {
        return d->metadata->assertContainsProperty(metakey);
    }
}

/** Return the list of metadata keys for the property with key 'key'

    \throw SireBase::missing_property
*/
QStringList Properties::metadataKeys(const PropertyName &key) const
{
    return this->allMetadata(key).propertyKeys();
}

/** Assert that this set contains the metadata property for the
    property 'key' with metakey 'metakey'

    \throw SireBase::missing_property
*/
void Properties::assertContainsMetadata(const PropertyName &key, const PropertyName &metakey) const
{
    this->assertContainsProperty(key);

    if (metakey.hasSource() and not this->allMetadata(key).hasProperty(metakey))
    {
        throw SireBase::missing_property(
            QObject::tr("There is no metadata with metakey \"%1\" for the property \"%2\". "
                        "Available metakeys are ( %3 ).")
                .arg(metakey.source(), key.toString(), metadataKeys(key).join(", ")),
            CODELOC);
    }
}

/** Return whether or not this contains a property with key 'key' */
bool Properties::hasProperty(const PropertyName &key) const
{
    if (this->isEmpty())
        return false;
    else if (key.hasValue() or d->properties.contains(key.source()) or key.hasDefaultValue())
        return true;
    else
    {
        auto linked = this->_getLink(key.source());
        return (not linked.isEmpty()) and d->properties.contains(linked);
    }
}

/** Return whether or not this contains the metadata with metakey 'metakey' */
bool Properties::hasMetadata(const PropertyName &metakey) const
{
    if (this->isEmpty())
        return false;
    else if (d->metadata.get() == 0)
        return false;
    else
        return d->metadata->hasProperty(metakey);
}

/** Return the property with key 'key'

    \throw SireBase::missing_property
*/
const Property &Properties::operator[](const PropertyName &key) const
{
    if (key.hasValue())
        return key.value();
    else
    {
        if (this->isEmpty())
        {
            if (key.hasDefaultValue())
                return key.value();
            else
                throw SireBase::missing_property(QObject::tr("There is no property with name \"%1\" "
                                                             "as there are no properties in this set.")
                                                     .arg(key.source()),
                                                 CODELOC);
        }

        QHash<QString, PropertyPtr>::const_iterator it = d->properties.constFind(key.source());

        if (it == d->properties.constEnd())
        {
            if (key.hasDefaultValue())
                return key.value();
            else
            {
                auto linked = this->_getLink(key.source());

                if (not linked.isEmpty())
                    it = d->properties.constFind(linked);

                if (it == d->properties.constEnd())
                    throw SireBase::missing_property(QObject::tr("There is no property with name \"%1\". "
                                                                 "Available properties are [ %2 ].")
                                                         .arg(key.source(), this->propertyKeys().join(", ")),
                                                     CODELOC);
            }
        }

        return *it;
    }
}

static Properties null_properties;

/** Return all of the metadata associated with this properties object */
const Properties &Properties::allMetadata() const
{
    if (this->isEmpty())
        return null_properties;
    else if (d->metadata.get() == 0)
        return null_properties;
    else
        return *(d->metadata);
}

/** Return the metadata for the property with key 'key'

    \throw SireBase::missing_property
*/
const Properties &Properties::allMetadata(const PropertyName &key) const
{
    if (key.hasValue())
        return null_properties;

    this->assertContainsProperty(key);

    if (d->props_metadata.get() == 0)
    {
        return null_properties;
    }
    else
    {
        auto it = d->props_metadata->constFind(key.source());

        if (it == d->props_metadata->constEnd())
            return null_properties;
        else
            return it.value();
    }
}

/** Return whether or not the property with key 'key' contains some
    metadata with the metakey 'metakey'

    \throw SireBase::missing_property
*/
bool Properties::hasMetadata(const PropertyName &key, const PropertyName &metakey) const
{
    this->assertContainsProperty(key);

    return metakey.hasValue() or this->allMetadata(key).hasProperty(metakey);
}

/** Return the property with key 'key' - note that if
    'key' specifies a value rather than a source, then the
    value contained in the key is returned

    \throw SireBase::missing_property
*/
const Property &Properties::property(const PropertyName &key) const
{
    return this->operator[](key);
}

/** Return the property with key 'key' - note that if
    'key' specifies a value rather than a source, then the
    value contained in the key is returned. If no such source
    exists, and there is no value in the key, then
    'default_value' is returned */
const Property &Properties::property(const PropertyName &key, const Property &default_value) const
{
    if (key.hasValue())
    {
        return key.value();
    }
    else if (this->isEmpty())
    {
        return default_value;
    }
    else
    {
        auto it = d->properties.constFind(key.source());

        if (it != d->properties.constEnd())
            return it.value();
        else
        {
            auto linked = this->_getLink(key.source());

            if (not linked.isEmpty())
            {
                it = d->properties.constFind(linked);

                if (it != d->properties.constEnd())
                    return it.value();
            }

            return default_value;
        }
    }
}

/** Return the metadata at metakey 'metakey' - note that if 'metakey'
    specifies a value rather than a source, then the value contained
    in the metakey is returned.

    \throw SireBase::missing_property
*/
const Property &Properties::metadata(const PropertyName &metakey) const
{
    if (this->isEmpty())
    {
        if (metakey.hasValue())
            return metakey.value();
    }
    else if (d->metadata.get() == 0)
    {
        if (metakey.hasValue())
            return metakey.value();
    }
    else
    {
        return d->metadata->property(metakey);
    }

    this->assertContainsMetadata(metakey);
    return metakey.value();
}

/** Return the metadata at metakey 'metakey' - note that if 'metakey'
    specifies a value rather than a source, then the value contained
    in the metakey is returned. If there is no such metadata, and no
    value is contained in the metakey, then 'default_value' is
    returned */
const Property &Properties::metadata(const PropertyName &metakey, const Property &default_value) const
{
    if (this->isEmpty())
    {
        if (metakey.hasValue())
            return metakey.value();
    }
    else if (d->metadata.get() == 0)
    {
        if (metakey.hasValue())
            return metakey.value();
    }
    else
    {
        return d->metadata->property(metakey, default_value);
    }

    return default_value;
}

/** Return the metadata at metakey 'metakey' that is associated with
    the property at key 'key'.

    \throw SireBase::missing_property
*/
const Property &Properties::metadata(const PropertyName &key, const PropertyName &metakey) const
{
    return this->allMetadata(key).property(metakey);
}

/** Return the metadata at metakey 'metakey' that is associated with
    the property at key 'key', or 'default_value' if there is no
    such metadata.
*/
const Property &Properties::metadata(const PropertyName &key, const PropertyName &metakey,
                                     const Property &default_value) const
{
    return this->allMetadata(key).property(metakey, default_value);
}

/** Set the property at key 'key' to have the value 'value'. This
    replaces any existing property at this key, and removes any
    existing metadata is 'clear_metadata' is true */
void Properties::setProperty(const QString &key, const Property &value, bool clear_metadata)
{
    if (key.isEmpty())
        throw SireError::invalid_arg(QObject::tr("You cannot insert a property with an empty key!"), CODELOC);

    if (this->isEmpty())
    {
        d = new PropertiesData();
    }

    d->properties.insert(key, value);

    // remove any link associated with this key name
    if (d->prop_links.get() != 0)
    {
        d->prop_links->remove(key);

        if (d->prop_links->isEmpty())
            d->prop_links.reset();
    }

    if (d->props_metadata.get() != 0)
    {
        if (clear_metadata or not d->props_metadata->contains(key))
        {
            d->props_metadata->remove(key);

            if (d->props_metadata->isEmpty())
                d->props_metadata.reset();
        }
    }
}

void Properties::setProperty(const QString &key, const Property &value)
{
    this->setProperty(key, value, false);
}

/** Set the metadata at metakey 'metakey' to have the value 'value'.
    This replaces any existing metadata with this metakey */
void Properties::setMetadata(const QString &metakey, const Property &value)
{
    if (this->isEmpty())
    {
        d = new PropertiesData();
    }

    if (d->metadata.get() == 0)
    {
        d->metadata.reset(new Properties());
    }

    d->metadata->setProperty(metakey, value);
}

/** Set the metadata at metakey 'metakey' for the property at key 'key'.
    This replaces any existing metadata for this key/metakey pair */
void Properties::setMetadata(const QString &key, const QString &metakey, const Property &value)
{
    this->assertContainsProperty(key);

    if (d->props_metadata.get() == 0)
    {
        d->props_metadata.reset(new QHash<QString, Properties>());
    }

    d->props_metadata->find(key)->setProperty(metakey, value);
}

/** Remove the property with key 'key' and all of its metadata */
void Properties::removeProperty(const QString &key)
{
    if (this->hasProperty(key))
    {
        d->properties.remove(key);

        if (d->props_metadata.get() != 0)
        {
            d->props_metadata->remove(key);

            if (d->props_metadata->isEmpty())
            {
                d->props_metadata.reset();
            }
        }

        if (d->prop_links.get() != 0)
        {
            QMutableHashIterator<QString, QString> it(*(d->prop_links));

            while (it.hasNext())
            {
                if (it.value() == key)
                    it.remove();
                else
                    it.next();
            }

            if (d->prop_links->isEmpty())
                d->prop_links.reset();
        }

        if (d->properties.isEmpty() and d->metadata.get() == 0)
        {
            d = 0;
        }
    }
}

/** Remove the metadata at metakey 'metakey' */
void Properties::removeMetadata(const QString &metakey)
{
    if (this->hasMetadata(metakey))
    {
        d->metadata->removeProperty(metakey);

        if (d->metadata->isEmpty())
        {
            d->metadata.reset();

            if (d->properties.isEmpty())
            {
                d = 0;
            }
        }
    }
}

/** Remove all of the top-level metadata */
void Properties::removeAllMetadata()
{
    if (this->isEmpty())
        return;

    if (d->metadata.get() != 0)
    {
        d->metadata->clear();

        if (d->properties.isEmpty())
            d = 0;
    }
}

/** Remove the metadata at metakey 'metakey' for the
    property at key 'key' */
void Properties::removeMetadata(const QString &key, const QString &metakey)
{
    if (this->hasMetadata(key, metakey))
    {
        d->props_metadata->find(key)->removeProperty(metakey);

        if (d->props_metadata->isEmpty())
        {
            d->props_metadata.reset();
        }
    }
}

/** Remove all of the metadata associated with the property at
    key 'key' */
void Properties::removeAllMetadata(const QString &key)
{
    if (this->hasProperty(key))
    {
        if (d->props_metadata.get() != 0)
        {
            d->props_metadata->remove(key);

            if (d->props_metadata->isEmpty())
                d->props_metadata.reset();
        }
    }
}

/** Completely clear this object of all properties and metadata */
void Properties::clear()
{
    d = 0;
}

/** Return the type name of the property at key 'key'

    \throw SireBase::missing_property
*/
const char *Properties::propertyType(const PropertyName &key) const
{
    return this->property(key).what();
}

/** Return the type name of the metadata at metakey 'metakey'

    \throw SireBase::missing_property
*/
const char *Properties::metadataType(const PropertyName &metakey) const
{
    return this->metadata(metakey).what();
}

/** Return the type name of the metadata at metakey 'metakey'
    for the property at key 'key'

    \throw SireBase::missing_property
*/
const char *Properties::metadataType(const PropertyName &key, const PropertyName &metakey) const
{
    return this->metadata(key, metakey).what();
}

/** Return the number of properties in this set */
int Properties::count() const
{
    if (this->isEmpty())
        return 0;
    else
        return d->properties.count();
}

/** Return the number of properties in this set */
int Properties::size() const
{
    return this->count();
}

/** Return the number of properties in this set */
int Properties::nProperties() const
{
    return this->count();
}

/** Return a string representation of this set of properties */
QString Properties::toString() const
{
    if (this->isEmpty())
        return QObject::tr("Properties::null");

    QStringList props;

    for (Properties::const_iterator it = this->constBegin(); it != this->constEnd(); ++it)
    {
        props.append(QObject::tr("    %1 => %2").arg(it.key(), it.value().read().toString()));
    }

    if (props.count() == 1)
        return QObject::tr("Properties( %1 )").arg(props.at(0).simplified());
    else
        return QObject::tr("Properties(\n%1\n)").arg(props.join(",\n"));
}

/** Internal function to return an editable property from the passed value.
 *  This returns 0 if the property doesn't exist
 */
Property *Properties::getEditableProperty(const QString &key)
{
    if (this->isEmpty())
        return 0;

    auto it = d->properties.find(key);

    if (it != d->properties.end())
        return it.value().data();
    else
    {
        auto linked = this->_getLink(key);

        if (not linked.isEmpty())
        {
            it = d->properties.find(linked);

            if (it != d->properties.end())
                return it.value().data();
        }

        return 0;
    }
}

/** Update the passed property to have the value 'value'. This does
 *  an in-place update on the existing property (which must have
 *  a compatible type). If 'auto-add' is true, then this will add
 *  the property if it doesn't exist. This returns whether or not
 *  a property was updated (or added)
 */
bool Properties::updateProperty(const QString &key, const Property &value,
                                bool auto_add)
{
    auto prop = this->getEditableProperty(key);

    if (prop == 0)
    {
        if (auto_add)
        {
            this->setProperty(key, value);
            return true;
        }
        else
            return false;
    }
    else
    {
        prop->copy(value);
        return true;
    }
}
