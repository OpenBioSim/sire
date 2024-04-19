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

#ifndef SIREBASE_PAGECACHE_H
#define SIREBASE_PAGECACHE_H

#include "sireglobal.h"

#include <QByteArray>

#include <memory>

SIRE_BEGIN_HEADER

namespace SireBase
{
    namespace detail
    {
        class CacheData;
        class PageData;
        class HandleData;
    }

    /** This class manages a swap cache of binary data that can be
     *  paged to and from disk. The cache can receive binary data
     *  of any size, and will automatically manage the paging of
     *  that data to and from disk as it is accessed.
     *
     *  You can create different caches, and have control over the maximum
     *  size of each cache page.
     *
     *  Note that deleting the cache will delete all data contained
     *  therein - including data paged to disk
     */
    class SIREBASE_EXPORT PageCache
    {
    public:
        PageCache(unsigned int page_size = 32 * 1024 * 1024);
        PageCache(const QString &cache_dir,
                  unsigned int page_size = 32 * 1024 * 1024);
        PageCache(std::shared_ptr<detail::CacheData> data);
        PageCache(const PageCache &other);
        ~PageCache();

        PageCache &operator=(const PageCache &other);

        const char *what() const;
        static const char *typeName();

        PageCache *clone() const;

        QString toString() const;

        QString cacheDir() const;

        unsigned int pageSize() const;
        unsigned int nPages() const;
        unsigned int nBytes() const;

        unsigned int size() const;

        bool isValid() const;
        bool isNull() const;

        void assertValid() const;

        static QString getStatistics();

        /** This is a page in the cache. This can hold multiple
         *  objects - the whole page is either resident in memory
         *  or cached to disk.
         */
        class SIREBASE_EXPORT Page
        {
            friend class detail::HandleData;

        public:
            Page();
            Page(std::shared_ptr<detail::PageData> data);
            Page(const Page &other);
            ~Page();

            Page &operator=(const Page &other);

            const char *what() const;
            static const char *typeName();

            Page *clone() const;

            QString toString() const;

            bool isValid() const;
            bool isNull() const;

            void assertValid() const;

            bool isResident() const;
            bool isCached() const;

            unsigned int nBytes() const;
            unsigned int size() const;

            unsigned int maxBytes() const;

            PageCache parent() const;

        protected:
            QByteArray fetch(int offset, int n_bytes) const;

        private:
            std::shared_ptr<detail::PageData> p;
        };

        /** This is a handle to a piece of data that has
         *  been added to the cache. This will either contain
         *  the actual data, or will hold the information
         *  necessary to retrieve that data from disk.
         *
         *  Data is removed from the cache when all handles
         *  to it are deleted
         */
        class SIREBASE_EXPORT Handle
        {
        public:
            Handle();
            Handle(std::shared_ptr<detail::HandleData> data);
            Handle(const Handle &other);
            ~Handle();

            Handle &operator=(const Handle &other);

            const char *what() const;
            static const char *typeName();

            Handle *clone() const;

            QString toString() const;

            Page page() const;
            QByteArray fetch() const;

            PageCache parent() const;

            bool isValid() const;
            bool isNull() const;

            void assertValid() const;

            unsigned int size() const;
            unsigned int nBytes() const;

            void clear();
            void reset();

        private:
            std::shared_ptr<detail::HandleData> h;
        };

        Handle store(const QByteArray &data);

    private:
        std::shared_ptr<detail::CacheData> d;
    };
}

SIRE_EXPOSE_CLASS(SireBase::PageCache)
SIRE_EXPOSE_CLASS(SireBase::PageCache::Page)
SIRE_EXPOSE_CLASS(SireBase::PageCache::Handle)

SIRE_END_HEADER

#endif
