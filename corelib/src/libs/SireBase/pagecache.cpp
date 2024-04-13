/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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
  *  You can contact the authors via the website
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#include "pagecache.h"

#include "SireError/errors.h"

#include <QDir>

#include <boost/noncopyable.hpp>

namespace SireBase
{
    namespace detail
    {
        class CacheData : public boost::noncopyable
        {
        public:
            CacheData(QString cachedir, int page_size);
            ~CacheData();

            QString cacheDir() const;

            PageCache::Handle cache(const QByteArray &data);

            int pageSize() const;

            int nPages() const;
            int nBytes() const;

        private:
            QString cache_dir;
            int page_size;
        };

        class PageData : public boost::noncopyable
        {
        public:
            PageData(int max_size);
            ~PageData();

            int maxBytes() const;
            int nBytes() const;
            int bytesRemaining() const;

            bool isResident() const;
            bool isCached() const;

            PageCache::Handle cache(const QByteArray &data);

            QByteArray fetch(int offset, int n_bytes) const;

            PageCache parent() const;

        private:
            std::weak_ptr<CacheData> c;
            int max_bytes;
            int nbytes;
        };

        class HandleData : public boost::noncopyable
        {
        public:
            HandleData();
            ~HandleData();

            QByteArray fetch() const;

            PageCache::Page page() const;

            int nBytes() const;

            bool isValid() const;
            bool isNull() const;

            PageCache parent() const;

        private:
            /** The page containing this data */
            PageCache::Page p;

            /** The offset of this data within the page */
            int offset;

            /** The number of bytes of data */
            int nbytes;
        };
    }
}

using namespace SireBase;
using namespace SireBase::detail;

///////
/////// Implementation of detail::CacheData
///////

CacheData::CacheData(QString c, int p)
    : cache_dir(c), page_size(p)
{
    if (c.simplified().isEmpty())
    {
        // by default, go into the current directory
        cache_dir = ".";
    }

    if (page_size < 1024)
    {
        throw SireError::invalid_arg(
            QObject::tr("Page size must be greater than 1024!"), CODELOC);
    }
}

CacheData::~CacheData()
{
}

PageCache::Handle CacheData::cache(const QByteArray &data)
{
    return PageCache::Handle();
}

QString CacheData::cacheDir() const
{
    return cache_dir;
}

int CacheData::pageSize() const
{
    return page_size;
}

int CacheData::nPages() const
{
    return 0;
}

int CacheData::nBytes() const
{
    return 0;
}

///////
/////// Implementation of detail::PageData
///////

PageData::PageData(int max_size)
    : max_bytes(max_size), nbytes(0)
{
    if (max_bytes < 1024)
    {
        throw SireError::invalid_arg(
            QObject::tr("Page size must be greater than 1024!"), CODELOC);
    }
}

PageData::~PageData()
{
}

int PageData::maxBytes() const
{
    return max_bytes;
}

int PageData::nBytes() const
{
    return nbytes;
}

int PageData::bytesRemaining() const
{
    return max_bytes - nbytes;
}

bool PageData::isResident() const
{
    return true;
}

bool PageData::isCached() const
{
    return false;
}

PageCache::Handle PageData::cache(const QByteArray &data)
{
    return PageCache::Handle();
}

QByteArray PageData::fetch(int offset, int n_bytes) const
{
    return QByteArray();
}

PageCache PageData::parent() const
{
    return PageCache(c.lock());
}

///////
/////// Implementation of detail::HandleData
///////

HandleData::HandleData()
    : offset(0), nbytes(0)
{
}

HandleData::~HandleData()
{
}

QByteArray HandleData::fetch() const
{
    return p.fetch(offset, nbytes);
}

PageCache::Page HandleData::page() const
{
    return p;
}

int HandleData::nBytes() const
{
    return nbytes;
}

bool HandleData::isValid() const
{
    return nbytes > 0;
}

bool HandleData::isNull() const
{
    return nbytes == 0;
}

PageCache HandleData::parent() const
{
    return p.parent();
}

///////
/////// Implementation of PageCache::Page
///////

PageCache::Page::Page()
    : p(nullptr)
{
}

PageCache::Page::Page(std::shared_ptr<detail::PageData> data)
    : p(data)
{
}

PageCache::Page::Page(const PageCache::Page &other)
    : p(other.p)
{
}

PageCache::Page::~Page()
{
}

PageCache::Page &PageCache::Page::operator=(const PageCache::Page &other)
{
    p = other.p;
    return *this;
}

const char *PageCache::Page::typeName()
{
    return "SireBase::PageCache::Page";
}

const char *PageCache::Page::what() const
{
    return PageCache::Page::typeName();
}

QString PageCache::Page::toString() const
{
    return QString("PageCache::Page");
}

void PageCache::Page::assertValid() const
{
    if (p == nullptr)
    {
        throw SireError::invalid_state(
            QObject::tr("Page object is null"), CODELOC);
    }
}

bool PageCache::Page::isValid() const
{
    return p != nullptr;
}

bool PageCache::Page::isNull() const
{
    return p == nullptr;
}

bool PageCache::Page::isResident() const
{
    assertValid();
    return p->isResident();
}

bool PageCache::Page::isCached() const
{
    assertValid();
    return p->isCached();
}

int PageCache::Page::nBytes() const
{
    assertValid();
    return p->nBytes();
}

int PageCache::Page::size() const
{
    return this->nBytes();
}

PageCache PageCache::Page::parent() const
{
    assertValid();
    return p->parent();
}

QByteArray PageCache::Page::fetch(int offset, int n_bytes) const
{
    assertValid();
    return p->fetch(offset, n_bytes);
}

///////
/////// Implementation of PageCache::Handle
///////

PageCache::Handle::Handle()
    : h(nullptr)
{
}

PageCache::Handle::Handle(std::shared_ptr<detail::HandleData> data)
    : h(data)
{
}

PageCache::Handle::Handle(const PageCache::Handle &other)
    : h(other.h)
{
}

PageCache::Handle::~Handle()
{
}

PageCache::Handle &PageCache::Handle::operator=(const PageCache::Handle &other)
{
    h = other.h;
    return *this;
}

const char *PageCache::Handle::typeName()
{
    return "SireBase::PageCache::Handle";
}

const char *PageCache::Handle::what() const
{
    return PageCache::Handle::typeName();
}

QString PageCache::Handle::toString() const
{
    return QString("PageCache::Handle");
}

void PageCache::Handle::assertValid() const
{
    if (h == nullptr)
    {
        throw SireError::invalid_state(
            QObject::tr("Handle object is null"), CODELOC);
    }
}

PageCache::Page PageCache::Handle::page() const
{
    assertValid();
    return Page(h->page());
}

QByteArray PageCache::Handle::fetch() const
{
    assertValid();
    return h->fetch();
}

PageCache PageCache::Handle::parent() const
{
    assertValid();
    return h->parent();
}

bool PageCache::Handle::isValid() const
{
    return h != nullptr;
}

bool PageCache::Handle::isNull() const
{
    return h == nullptr;
}

int PageCache::Handle::size() const
{
    return this->nBytes();
}

int PageCache::Handle::nBytes() const
{
    assertValid();
    return h->nBytes();
}

void PageCache::Handle::clear()
{
    h.reset();
}

void PageCache::Handle::reset()
{
    h.reset();
}

///////
/////// Implementation of PageCache
///////

PageCache::PageCache(int page_size)
    : d(new CacheData(QString("."), page_size))
{
}

PageCache::PageCache(const QString &cachedir, int page_size)
    : d(new CacheData(cachedir, page_size))
{
}

PageCache::PageCache(std::shared_ptr<detail::CacheData> data)
    : d(data)
{
}

PageCache::PageCache(const PageCache &other)
    : d(other.d)
{
}

PageCache::~PageCache()
{
}

PageCache &PageCache::operator=(const PageCache &other)
{
    d = other.d;
    return *this;
}

const char *PageCache::typeName()
{
    return "SireBase::PageCache";
}

const char *PageCache::what() const
{
    return PageCache::typeName();
}

QString PageCache::toString() const
{
    return QString("PageCache");
}

void PageCache::assertValid() const
{
    if (d == nullptr)
    {
        throw SireError::invalid_state(
            QObject::tr("PageCache object is null"), CODELOC);
    }
}

QString PageCache::cacheDir() const
{
    assertValid();
    return d->cacheDir();
}

int PageCache::pageSize() const
{
    assertValid();
    return d->pageSize();
}

int PageCache::nPages() const
{
    assertValid();
    return d->nPages();
}

int PageCache::nBytes() const
{
    assertValid();
    return d->nBytes();
}

int PageCache::size() const
{
    return this->nBytes();
}

bool PageCache::isValid() const
{
    return d != nullptr;
}

bool PageCache::isNull() const
{
    return d == nullptr;
}

PageCache::Handle PageCache::cache(const QByteArray &data)
{
    assertValid();
    return Handle(d->cache(data));
}
