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

#include "SireBase/parallel.h"
#include "SireBase/console.h"
#include "SireBase/atexit.h"

#include <QDir>
#include <QThread>
#include <QQueue>
#include <QHash>
#include <QMutex>
#include <QTemporaryFile>
#include <QTemporaryDir>
#include <QAtomicInt>

#include <boost/noncopyable.hpp>

#include <cstring> // needed for memcpy

namespace SireBase
{
    namespace detail
    {
        /** This class holds a page of data that has been restored from disk */
        class RestoredPage : public boost::noncopyable
        {
        public:
            RestoredPage(const QByteArray &data);
            ~RestoredPage();

            QByteArray fetch(unsigned int offset, unsigned int n_bytes) const;

            int timeToLive() const;

        private:
            /** The restored data */
            QByteArray d;

            /** The time to live - this is reset every time the
             *  data is accessed via the fetch method */
            QAtomicInt ttl;
        };

        /** This class handles the movement of pages between disk
         *  and memory
         */
        class PageHandler : boost::noncopyable
        {
        public:
            PageHandler(QString cache_dir_template);
            ~PageHandler();

            QPair<QTemporaryFile *, std::shared_ptr<RestoredPage>> store(const QByteArray &data);
            std::shared_ptr<RestoredPage> restore(QTemporaryFile &pagefile);

            QString path() const;

            void cleanUpOnExit();

            static void setMaxResidentPages(unsigned int n_pages);
            static unsigned int maxResidentPages();

        private:
            void _lkr_addToRestored(const std::shared_ptr<RestoredPage> &page);

            /** The maximum number of pages left resident in memory
             *  per cache
             */
            static unsigned int max_resident_pages;

            /** Mutex to protect access to the data of this class */
            QMutex mutex;

            /** List of all pages restored by this handler */
            QList<std::shared_ptr<RestoredPage>> restored_pages;

            /** Checksum of the data in all temporary page files */
            QHash<QTemporaryFile *, quint16> checksums;

            /** Pointer to the temporary directory for all pagefiles */
            QTemporaryDir *cache_dir;
        };

        /** This class holds all of the data for a PageCache, providing
         *  a cache where data is stored, returning a handle to the
         *  stored data, through which it can be fetched. The data
         *  is internally divided into pages, which are pushed to disk
         *  so that memory usage is kept to a minimum.
         */
        class CacheData : public QThread
        {
        public:
            CacheData(QString cachedir, unsigned int page_size);
            ~CacheData();

            static QString getStatistics();

            static void setMaxPageSize(unsigned int size, bool update_existing);
            static unsigned int maxPageSize();

            QString cacheDir() const;

            PageCache::Handle store(const QByteArray &data);

            unsigned int pageSize() const;

            unsigned int nPages() const;
            unsigned int nBytes() const;

            void setSelf(const std::shared_ptr<CacheData> &cache);

            static void cleanUpOnExit();

        protected:
            void run();

        private:
            void enqueue(const std::shared_ptr<HandleData> &handle);
            std::shared_ptr<PageHandler> getPageHandler();

            /** Mutex protecting the list of all caches */
            static QMutex caches_mutex;

            /** The list of all caches - weak pointers so that they
             *  can be deleted when the last reference is removed
             */
            static QList<std::weak_ptr<CacheData>> caches;

            /** The current default maximum page size */
            static unsigned int max_page_size;

            /** Mutex protected the queue of Handles containing data
             *  that should be pushed to the cache
             */
            QMutex queue_mutex;

            /** Queue of handles that are yet to be cached */
            QQueue<std::weak_ptr<HandleData>> queue;

            /** Mutex to protect access to the list of pages
             *  associated with this cache
             */
            QMutex page_mutex;

            /** The current page of data being filled - weak pointer
             *  so that it is deleted if no longer needed
             */
            std::weak_ptr<PageData> current_page;

            /** The list of filled pages of data, which have been frozen
             *  to disk. Weak pointers so that they are deleted if no longer
             *  needed
             */
            QList<std::weak_ptr<PageData>> pages;

            /** Weak point to self */
            std::weak_ptr<CacheData> self;

            /** Weak pointer to the page handler for this cache.
             *  This is passed to pages when they are frozen so that
             *  the handler can be used to restore the data.
             *  Weak pointer to that is is automatically deleted when
             *  all pages are deleted
             */
            std::weak_ptr<PageHandler> page_handler;

            /** Template used for the cache directory - this is in a format
             *  that is recognised by QTemporaryDir
             */
            QString cache_dir_template;

            /** The default maximum page size for pages. This is a guide,
             *  i.e. items larger than this size will be cached into their
             *  own page. Also, a new page will be automatically created
             *  if there is insufficient space on the current page
             */
            int page_size;

            /** Flag set when the underlying thread is exiting. This is used
             *  to prevent race conditions
             */
            bool exiting;
        };

        /** This class holds the implementation of a single page of data */
        class PageData : public boost::noncopyable
        {
        public:
            PageData(unsigned int max_size, const std::shared_ptr<CacheData> &cache);
            ~PageData();

            unsigned int maxBytes() const;
            unsigned int nBytes() const;
            unsigned int bytesRemaining() const;

            bool isResident() const;
            bool isCached() const;

            unsigned int store(const QByteArray &data);

            void freeze(std::shared_ptr<PageHandler> handler);

            QByteArray fetch(unsigned int offset, unsigned int n_bytes) const;

            PageCache parent() const;

            void cleanUpOnExit();

        private:
            /** Weak pointer to the parent cache - weak so that the
             *  cache is automatically deleted if the user releases
             *  all references to it
             */
            std::weak_ptr<CacheData> c;

            /** Mutex to protect access to the data of this page */
            QMutex mutex;

            /** The handler used to restore this page from disk.
             *  This is null if the page has not yet been frozen
             */
            std::shared_ptr<PageHandler> page_handler;

            /** Weak pointer to the restored page. This is weak so
             *  that, if the page handler decides to remove the page
             *  from memory, then this pointer will be automatically
             *  set to null in a thread-safe way
             */
            std::weak_ptr<RestoredPage> restored_page;

            /** The maximum number of bytes allowed in the page.
             *  This is set equal to the current number of bytes
             *  when the page is frozen
             */
            unsigned int max_bytes;

            /** The current number of bytes in the page */
            unsigned int nbytes;

            /** Pointer to the temporary file used to hold the frozen
             *  page
             */
            QTemporaryFile *cache_file;

            /** Pointer to the page of data - this is only used when the
             *  page is being filled. It is set to null when the page
             *  is frozen to disk
             */
            char *d;
        };

        class HandleData : public boost::noncopyable
        {
        public:
            enum State
            {
                EMPTY = 0,
                DATA_IN_HANDLE = 1,
                DATA_ON_PAGE = 2
            };

            HandleData(const QByteArray &data);
            ~HandleData();

            QByteArray fetch() const;

            PageCache::Page page() const;

            unsigned int nBytes() const;

            bool isValid() const;
            bool isNull() const;

            PageCache parent() const;

            void setPage(const PageCache::Page &page, int offset);

        private:
            /** The page containing this data */
            PageCache::Page p;

            /** Mutex to lock access to this data */
            QMutex mutex;

            /** The data being cached - this is stored here until it is
                property added to the cache */
            QByteArray d;

            /** The current state of the handle */
            State state;

            /** The offset of this data within the page */
            unsigned int offset;

            /** The number of bytes of data */
            unsigned int nbytes;
        };
    }
}

using namespace SireBase;
using namespace SireBase::detail;

///////
/////// Implementation of detail::RestoredPage
///////

/** Construct to hold the restored data in 'data'. By default,
 *  the ttl starts at 100
 */
RestoredPage::RestoredPage(const QByteArray &data)
    : d(data), ttl(100)
{
}

/** Destructor */
RestoredPage::~RestoredPage()
{
}

/** Return the time to live for this page, and decrement it by 1 */
int RestoredPage::timeToLive() const
{
    // this uses an atomic integer, so is thread-safe
    return const_cast<RestoredPage *>(this)->ttl.fetchAndAddRelaxed(-1);
}

/** Fetch 'n_bytes' bytes of data starting at 'offset' */
QByteArray RestoredPage::fetch(unsigned int offset, unsigned int n_bytes) const
{
    if (offset + n_bytes > static_cast<unsigned int>(d.size()))
    {
        throw SireError::invalid_arg(
            QObject::tr("Impossible to fetch %1 bytes starting at "
                        "offset %2 from a page of only %3 bytes")
                .arg(n_bytes)
                .arg(offset)
                .arg(d.size()),
            CODELOC);
    }

    // reset the TTL to 100 - this uses an atomic integer, so is thread-safe
    const_cast<RestoredPage *>(this)->ttl.storeRelaxed(100);

    // return the requested data
    return QByteArray(d.constData() + offset, n_bytes);
}

///////
/////// Implementation of detail::PageHandler
///////

/** Construct a PageHandler that will handle freezing and restoring
 *  pages from a temporary directory following the passed template
 */
PageHandler::PageHandler(QString cache_dir_template)
    : cache_dir(new QTemporaryDir(cache_dir_template))
{
    if (not cache_dir->isValid())
    {
        auto message = QObject::tr("Failed to create cache directory %1. %2")
                           .arg(cache_dir_template)
                           .arg(cache_dir->errorString());

        Console::error(message);
        throw SireError::io_error(message, CODELOC);
    }
}

/** Destructor */
PageHandler::~PageHandler()
{
    delete cache_dir;
}

/** Clean up the PageHandler (called on exit) */
void PageHandler::cleanUpOnExit()
{
    QMutexLocker lkr(&mutex);
    delete cache_dir;
    cache_dir = 0;
    checksums.clear();
    restored_pages.clear();
}

/** Return the path to the cache directory */
QString PageHandler::path() const
{
    return cache_dir->path();
}

// this would be a maximum of 256 MB resident per cache
unsigned int PageHandler::max_resident_pages = 32;

/** Set the maximum number of pages that can be resident in memory
 *  at any one time
 */
void PageHandler::setMaxResidentPages(unsigned int n_pages)
{
    if (n_pages < 1)
    {
        Console::warning(QObject::tr("Setting maximum resident pages to the minimum of 1."));
        n_pages = 1;
    }
    else if (n_pages > 256)
    {
        Console::warning(QObject::tr("Setting maximum resident pages to the maximum of 256."));
        n_pages = 256;
    }

    max_resident_pages = n_pages;
}

/** Return the current maximum number of pages that can be resident
 *  in memory at any one time
 */
unsigned int PageHandler::maxResidentPages()
{
    return max_resident_pages;
}

/** Add the passed restored page to the cache */
void PageHandler::_lkr_addToRestored(const std::shared_ptr<RestoredPage> &restored)
{
    // check to see if we need to replace an old page
    auto max_resident = max_resident_pages;

    if (static_cast<unsigned int>(restored_pages.count()) < max_resident)
    {
        restored_pages.append(restored);
        return;
    }

    // we need to replace a page...
    int smallest_index = -1;
    int smallest_value = -1;

    for (int i = 0; i < restored_pages.size(); i++)
    {
        if (restored_pages[i].get() == nullptr)
        {
            // this page has already been deleted, so should
            // be removed first
            smallest_index = i;
        }
        else
        {
            // find the oldest page to restore - do this by reducing
            // the TTL for each page each time we go through this loop,
            // and then find the first page with the smallest TTL
            int val = restored_pages[i]->timeToLive();

            if (smallest_index == -1 or smallest_value > val)
            {
                smallest_index = i;
                smallest_value = val;
            }
        }
    }

    restored_pages[smallest_index] = restored;
}

/** Store the passed page of data to a temporary page cache file.
 *  This returns a pointer to the created temporary file,
 *  which is owned by the caller (in this case, the PageData)
 */
QPair<QTemporaryFile *, std::shared_ptr<RestoredPage>> PageHandler::store(const QByteArray &data)
{
    if (cache_dir == 0)
    {
        // we are likely being shut down, so return a null pointer
        return QPair<QTemporaryFile *, std::shared_ptr<RestoredPage>>(0, 0);
    }

    QTemporaryFile *cache_file = new QTemporaryFile(cache_dir->filePath("page_XXXXXX.bin"));

    try
    {
        if (not cache_file->open())
        {
            Console::error(QObject::tr("Failed to open temporary file %1. %2")
                               .arg(cache_file->fileName())
                               .arg(cache_file->errorString()));
            throw SireError::file_error(*cache_file, CODELOC);
        }

        // create a handle to the restored page
        auto restored_page = std::make_shared<RestoredPage>(data);

        // compress the data and write it to a temporary file
        QByteArray compressed_data = qCompress(data, 9);

        auto checksum = qChecksum(compressed_data.constData(),
                                  compressed_data.size());

        if (compressed_data.size() != cache_file->write(compressed_data))
        {
            throw SireError::file_error(*cache_file, CODELOC);
        }

        cache_file->close();

        QMutexLocker lkr(&mutex);

        checksums.insert(cache_file, checksum);

        // save this restored page so that it is not immediately deleted
        this->_lkr_addToRestored(restored_page);

        return QPair<QTemporaryFile *, std::shared_ptr<RestoredPage>>(cache_file, restored_page);
    }
    catch (...)
    {
        delete cache_file;
        throw;
    }
}

/** Restore the QByteArray held in the passed pagefile, and
 *  add it to the list of restored pages. If there are more than
 *  5 restored pages, then the oldest page is replaced with the
 *  new page (oldest based on looking at the TTL values of each page)
 */
std::shared_ptr<RestoredPage> PageHandler::restore(QTemporaryFile &pagefile)
{
    if (not pagefile.open())
    {
        Console::error(QObject::tr("Failed to open temporary file %1. %2")
                           .arg(pagefile.fileName())
                           .arg(pagefile.errorString()));
        throw SireError::file_error(pagefile, CODELOC);
    }

    // read the raw data
    QByteArray compressed_data = pagefile.readAll();
    pagefile.close();

    // make sure that the checksum matches what we calculated when we
    // stored the data to disk
    QMutexLocker lkr(&mutex);

    auto checksum = checksums.value(&pagefile, 0);
    auto readsum = qChecksum(compressed_data.constData(), compressed_data.size());

    if (checksum != readsum)
    {
        QString message = QObject::tr("Checksum failed on page data for %1: %2 versus %3")
                              .arg(pagefile.fileName())
                              .arg(checksum)
                              .arg(readsum);

        Console::error(message);

        throw SireError::invalid_state(message, CODELOC);
    }

    // uncompress the data
    QByteArray decompressed_data = qUncompress(compressed_data);
    compressed_data = QByteArray();

    // create a new RestoredPage object for this data
    auto restored = std::make_shared<RestoredPage>(decompressed_data);
    decompressed_data = QByteArray();

    this->_lkr_addToRestored(restored);

    return restored;
}

///////
/////// Implementation of detail::CacheData
///////

QMutex CacheData::caches_mutex;
QList<std::weak_ptr<CacheData>> CacheData::caches;

/** Default to 8 MB pages - this seems to be a good balance for
 *  not having pages so large that they take a noticeable time to
 *  load from disk
 */
unsigned int CacheData::max_page_size = 8 * 1024 * 1024;

/** Clean-up function that will remove all temporary files etc when
 *  the program exits, and will interupt all of the background threads
 *  so that they will (hopefully) exit
 */
void CacheData::cleanUpOnExit()
{
    QMutexLocker lkr(&caches_mutex);

    for (auto &cache : caches)
    {
        auto c = cache.lock();

        if (c.get() != nullptr)
        {
            // interupt the thread so that it will stop
            c->requestInterruption();

            // clean up all of the pages associated with this cache
            QMutexLocker lkr2(&(c->page_mutex));

            for (auto &page : c->pages)
            {
                auto p = page.lock();

                if (p.get() != nullptr)
                {
                    p->cleanUpOnExit();
                }
            }

            lkr2.unlock();

            auto handler = c->page_handler.lock();

            if (handler.get() != nullptr)
            {
                handler->cleanUpOnExit();
            }
        }
    }

    // wait for all of the threads to stop, and if they don't
    // then kill them
    for (auto &cache : caches)
    {
        auto c = cache.lock();

        if (c.get() != nullptr)
        {
            if (QThread::currentThread() != c.get())
            {
                bool is_finished = false;

                // give the thread a maximum of 10 attempts to finish
                for (int i = 0; i < 10; i++)
                {
                    // the thread didn't finish
                    if (not c->wait(100))
                    {
                        qDebug() << "Waiting for cache thread to finish:" << c->cacheDir();
                    }
                    // the thread finished
                    else
                    {
                        is_finished = true;
                        break;
                    }
                }

                // the thread still didn't finish, so kill it
                if (not is_finished)
                {
                    qDebug() << "Terminating cache thread:" << c->cacheDir();
                    c->terminate();
                }
            }
        }
    }
}

/** This is the clean-up function that is registered with the
 *  atexit function. It will call the cleanUpOnExit function
 *  of CacheData
 */
void clean_up_pagecache()
{
    CacheData::cleanUpOnExit();
}

static const RegisterExitFunction clean_up_pagecache_instance(clean_up_pagecache);

/** Set the maximum page size to the passed value */
void CacheData::setMaxPageSize(unsigned int size, bool update_existing)
{
    if (size < 1024)
    {
        Console::warning(QObject::tr("Setting page size to the minimum of 1024 bytes."));
        size = 1024;
    }
    else if (size > 128 * 1024 * 1024)
    {
        Console::warning(QObject::tr("Setting page size to the maximum of 128 MB."));
        size = 128 * 1024 * 1024;
    }

    QMutexLocker lkr(&caches_mutex);
    max_page_size = size;

    if (update_existing)
    {
        for (auto &cache : caches)
        {
            auto c = cache.lock();

            if (c.get() != nullptr)
            {
                c->page_size = size;
            }
        }
    }
}

/** Return the current maximum page size */
unsigned int CacheData::maxPageSize()
{
    QMutexLocker lkr(&caches_mutex);
    return max_page_size;
}

/** Set the self pointer for this cache */
void CacheData::setSelf(const std::shared_ptr<CacheData> &cache)
{
    if (cache.get() == nullptr)
    {
        return;
    }

    if (self.lock().get() != nullptr)
    {
        QString message = QObject::tr("Cache already has a self pointer!");
        Console::error(message);
        throw SireError::invalid_state(message, CODELOC);
    }

    if (cache.get() != this)
    {
        QString message = QObject::tr("Cache self pointer does not match this cache!");
        Console::error(message);
        throw SireError::invalid_state(message, CODELOC);
    }

    QMutexLocker lkr(&caches_mutex);

    if (self.lock().get() != nullptr)
    {
        QString message = QObject::tr("Cache already has a self pointer!");
        Console::error(message);
        throw SireError::invalid_state(message, CODELOC);
    }

    caches.append(cache);
    self = cache;
}

/** Return a string giving all of the statistics for all caches */
QString CacheData::getStatistics()
{
    QMutexLocker lkr(&caches_mutex);
    QList<std::weak_ptr<CacheData>> local_caches = caches;
    lkr.unlock();

    QString stats;

    for (auto &cache : local_caches)
    {
        auto c = cache.lock();

        if (c.get() != nullptr)
        {
            stats += QString("Cache: %1\n")
                         .arg(c->cacheDir());

            QMutexLocker lkr2(&(c->page_mutex));

            auto current = c->current_page.lock();

            int total_bytes = 0;

            if (current.get() != nullptr)
            {
                stats += QString("  Current Page: %1 KB : is_resident %2\n")
                             .arg(current->nBytes() / 1024.0)
                             .arg(current->isResident());

                total_bytes += current->nBytes();
            }

            int page_count = 0;

            for (auto &page : c->pages)
            {
                auto p = page.lock();

                if (p.get() != nullptr)
                {
                    page_count += 1;

                    stats += QString("  Page %1: %2 KB : is_resident %3\n")
                                 .arg(page_count)
                                 .arg(p->nBytes() / 1024.0)
                                 .arg(p->isResident());

                    total_bytes += p->nBytes();
                }
            }

            stats += QString("  Total size: %1 KB\n")
                         .arg(total_bytes / 1024.0);
        }
    }

    return stats;
}

/** Construct the data for a single cache, using the passed
 *  template for the cache directory and the passed suggested
 *  maximum page size
 */
CacheData::CacheData(QString c, unsigned int p)
    : page_size(p), exiting(false)
{
    if (c.simplified().isEmpty())
    {
        // by default, go into the current directory
        c = QDir(PageCache::rootDirectory()).filePath("temp_XXXXXX");
    }

    if (page_size < 1024)
    {
        Console::warning(QObject::tr("Setting page size to the minimum of 1024 bytes."));
        page_size = 1024;
    }
    else if (page_size > 128 * 1024 * 1024)
    {
        Console::warning(QObject::tr("Setting page size to the maximum of 128 MB."));
        page_size = 128 * 1024 * 1024;
    }

    cache_dir_template = c;

    // make sure that we can actually create a cache directory in this
    // space - do this by creating a test directory and then deleting it
    QTemporaryDir test_dir(c);

    if (not test_dir.isValid())
    {
        throw SireError::io_error(QObject::tr(
                                      "Failed to create cache directory %1. %2")
                                      .arg(c)
                                      .arg(test_dir.errorString()),
                                  CODELOC);
    }
}

/** Destructor */
CacheData::~CacheData()
{
    // tell the thread to stop
    this->requestInterruption();

    if (QThread::currentThread() != this)
    {
        // wait for the thread to stop
        this->wait();
    }
    else
    {
        // we are in the thread, so we can't wait for it to stop
        // (use qWarning, as at this point the Console may be gone...)
        qWarning() << "CacheData is deleting itself!" << this->cacheDir();
    }
}

/** Return the PageHandler for this cache - this automatically creates
 *  a new one if it doesn't exist, or if the old one was lost
 *  (happens if all pages are deleted)
 *
 *  Note that this function is ONLY called from the background running
 *  thread, so does not need to be (and isn't) thread-safe
 */
std::shared_ptr<PageHandler> CacheData::getPageHandler()
{
    auto handler = page_handler.lock();

    if (handler.get() == nullptr)
    {
        handler = std::make_shared<PageHandler>(cache_dir_template);
        page_handler = handler;
    }

    return handler;
}

/** Store the data in the cache, returning a handle to the data */
PageCache::Handle CacheData::store(const QByteArray &data)
{
    auto handle = std::make_shared<HandleData>(data);
    this->enqueue(handle);
    return PageCache::Handle(handle);
}

/** Return the cache directory for this cache */
QString CacheData::cacheDir() const
{
    auto handler = page_handler.lock();

    if (handler.get() == nullptr)
    {
        return cache_dir_template;
    }
    else
    {
        return handler->path();
    }
}

/** Return the suggested maximum page size for this cache */
unsigned int CacheData::pageSize() const
{
    return page_size;
}

/** Return the number of pages in this cache */
unsigned int CacheData::nPages() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&page_mutex));
    unsigned int npages = 0;

    if (current_page.lock().get() != nullptr)
    {
        npages += 1;
    }

    for (auto &page : pages)
    {
        if (page.lock().get() != nullptr)
        {
            npages += 1;
        }
    }

    return npages;
}

/** Return the number of bytes in this cache */
unsigned int CacheData::nBytes() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&page_mutex));

    unsigned int nbytes = 0;

    auto current = current_page.lock();

    if (current.get() != nullptr)
    {
        nbytes += current->nBytes();
    }

    for (auto &page : pages)
    {
        auto p = page.lock();

        if (p.get() != nullptr)
        {
            nbytes += p->nBytes();
        }
    }

    return nbytes;
}

/** Internal function that is used to add a handler to the queue
 *  to be processed by the background thread and copied to a page
 */
void CacheData::enqueue(const std::shared_ptr<HandleData> &handle)
{
    if (handle.get() == nullptr)
    {
        return;
    }

    QMutexLocker lkr(&queue_mutex);

    if (this->exiting)
    {
        // this is the race condition where the current thread loop
        // is in the process of exiting, but we need to wait for that
        // to complete
        if (QThread::currentThread() != this)
            this->wait();

        this->exiting = false;
    }

    queue.enqueue(handle);

    // Start the background thread if it is not running
    // (it automatically stops if there is nothing to do)
    if (not this->isRunning())
    {
        this->start();
    }
}

/** The background thread that processes the queue of data to be
 *  added to the cache
 */
void CacheData::run()
{
    try
    {
        int empty_count = 0;

        while (true)
        {
            if (this->isInterruptionRequested())
            {
                // stop what we are doing
                break;
            }

            std::shared_ptr<HandleData> handle;
            bool have_item = false;

            // pull an item off the queue
            {
                QMutexLocker lkr(&queue_mutex);

                if (this->isInterruptionRequested())
                {
                    // stop what we are doing
                    break;
                }

                while (not queue.isEmpty())
                {
                    handle = queue.dequeue().lock();

                    if (handle.get() != nullptr)
                    {
                        have_item = true;
                        break;
                    }
                }
            }

            if (have_item)
            {
                // get the data
                QByteArray data = handle->fetch();

                const unsigned int n_bytes = static_cast<unsigned int>(data.size());

                if (n_bytes >= page_size)
                {
                    // this is bigger than a page, so needs to have its
                    // own page!
                    auto page = std::make_shared<PageData>(n_bytes, this->self.lock());
                    QMutexLocker lkr(&page_mutex);
                    this->pages.append(page);
                    lkr.unlock();
                    auto offset = page->store(data);

                    // straight write this to disk and prevent any changes
                    page->freeze(this->getPageHandler());

                    // return a handle to the new page
                    handle->setPage(PageCache::Page(page), offset);
                }
                else if (n_bytes != 0)
                {
                    // make sure we have a current page...
                    auto current = current_page.lock();

                    if (current.get() == nullptr)
                    {
                        current = std::make_shared<PageData>(page_size, this->self.lock());
                        QMutexLocker lkr(&page_mutex);
                        current_page = current;
                    }

                    // is there space left on the current page
                    if (current->bytesRemaining() >= n_bytes)
                    {
                        auto offset = current->store(data);
                        handle->setPage(PageCache::Page(current), offset);
                    }
                    else
                    {
                        // freeze this page to prevent any further changes
                        current->freeze(this->getPageHandler());

                        // add the current page to the list of old pages
                        QMutexLocker lkr(&page_mutex);
                        this->pages.append(current_page);

                        // we need to create a new page
                        current = std::make_shared<PageData>(page_size, this->self.lock());
                        current_page = current;
                        lkr.unlock();

                        auto offset = current->store(data);
                        handle->setPage(PageCache::Page(current), offset);
                    }
                }

                empty_count = 0;

                // we've finished with 'handle' - release the shared
                // pointer so we aren't holding onto it for longer than we need
                handle.reset();
            }

            if (this->isInterruptionRequested())
            {
                // stop what we are doing
                break;
            }

            bool is_empty = true;

            if (have_item)
            {
                // check to see if there is anything else to process
                QMutexLocker lkr(&queue_mutex);
                is_empty = queue.isEmpty();
            }

            if (is_empty)
            {
                if (this->isInterruptionRequested())
                {
                    // stop what we are doing
                    break;
                }

                empty_count += 1;

                if (empty_count > 10)
                {
                    // we have been idle for a while, so we can stop
                    this->exiting = true;
                    break;
                }

                // sleep for a bit
                this->msleep(100);
            }
        }
    }
    catch (const SireError::exception &e)
    {
        Console::error(e.error());
    }
    catch (const std::exception &e)
    {
        Console::error(e.what());
    }
    catch (...)
    {
        Console::error("Unknown error in CacheData::run()");
    }
}

///////
/////// Implementation of detail::PageData
///////

/** Construct an empty new page of data of specified maximum size,
 *  associated with the specified cache
 */
PageData::PageData(unsigned int max_size, const std::shared_ptr<CacheData> &cache)
    : c(cache), max_bytes(max_size), nbytes(0), cache_file(0), d(0)
{
    if (max_bytes < 1024)
    {
        Console::warning(QObject::tr("Setting page size to the minimum of 1024 bytes."));
        max_bytes = 1024;
    }
    else if (max_bytes > 128 * 1024 * 1024)
    {
        Console::warning(QObject::tr("Setting page size to the maximum of 128 MB."));
        max_bytes = 128 * 1024 * 1024;
    }

    d = new char[max_bytes];
}

/** Destructor */
PageData::~PageData()
{
    delete cache_file;
    delete[] d;
}

/** Called when the program is exiting - makes sure that page file is
 *  deleted
 */
void PageData::cleanUpOnExit()
{
    QMutexLocker lkr(&mutex);

    if (cache_file != 0)
    {
        cache_file->remove();
        delete cache_file;
        cache_file = 0;
    }

    if (d != 0)
    {
        delete[] d;
        d = 0;
    }
}

/** Return the maximum number of bytes that can be stored in this page */
unsigned int PageData::maxBytes() const
{
    if (d == 0)
    {
        return nbytes;
    }
    else
    {
        return max_bytes;
    }
}

/** Return the number of bytes currently stored in this page */
unsigned int PageData::nBytes() const
{
    return nbytes;
}

/** Return the number of bytes remaining in this page */
unsigned int PageData::bytesRemaining() const
{
    if (d == 0)
    {
        return 0;
    }
    else
    {
        return max_bytes - nbytes;
    }
}

/** Return true if the page is resident in memory */
bool PageData::isResident() const
{
    return d != 0 or restored_page.lock().get() != nullptr;
}

/** Return true if the page is cached to disk */
bool PageData::isCached() const
{
    return cache_file != 0;
}

/** Store the data in the page, returning the offset at which
 *  the data is stored. This function is only ever called by
 *  the background thread of a CacheData, so does not need to
 *  be (and isn't) thread-safe
 */
unsigned int PageData::store(const QByteArray &data)
{
    // this test will fail if the page is frozen
    if (static_cast<unsigned int>(data.size()) > this->bytesRemaining())
    {
        QString message = QObject::tr("Data is too large to fit on this page!");
        Console::error(message);
        throw SireError::invalid_arg(
            QObject::tr("Data is too large to fit on this page!"), CODELOC);
    }

    // only the thread calling this function will change d, so it is
    // save to assume that d will not change during this function

    if (d == 0)
    {
        QString message = QObject::tr("Page is already frozen!");
        Console::error(message);
        throw SireError::invalid_state(
            QObject::tr("Page is already frozen!"), CODELOC);
    }

    // copy the data into the page
    std::memcpy(d + nbytes, data.constData(), data.size());

    // update the number of bytes in the page and get the offset
    unsigned int offset = nbytes;
    nbytes += data.size();

    return offset;
}

/** Fetch 'n_bytes' bytes of data starting at 'offset' from this page */
QByteArray PageData::fetch(unsigned int offset, unsigned int n_bytes) const
{
    if (offset + n_bytes > nbytes)
    {
        throw SireError::invalid_arg(
            QObject::tr("Impossible to fetch %1 bytes starting at "
                        "offset %2 from a page of only %3 bytes")
                .arg(n_bytes)
                .arg(offset)
                .arg(nbytes),
            CODELOC);
    }

    auto restored = restored_page.lock();

    if (restored.get() != nullptr)
    {
        return restored->fetch(offset, n_bytes);
    }

    // need to lock because the background thread may be changing
    // the data (moving d to the page)
    QMutexLocker lkr(&(const_cast<PageData *>(this)->mutex));

    if (d == 0)
    {
        if (page_handler.get() == nullptr or cache_file == 0)
        {
            QString message = QObject::tr("Page has not been frozen to disk?");
            Console::error(message);
            return QByteArray();
        }

        auto restored = page_handler->restore(*cache_file);

        const_cast<PageData *>(this)->restored_page = restored;
        lkr.unlock();

        return restored->fetch(offset, n_bytes);
    }
    else
    {
        return QByteArray(d + offset, n_bytes);
    }
}

/** Freeze the page to disk, using the passed handler. This function
 *  will only ever be called by the background thread of a CacheData,
 *  so does not need to be (and isn't) thread-safe
 */
void PageData::freeze(std::shared_ptr<PageHandler> handler)
{
    // page is already frozen, or something else weird
    if (d == 0 or nbytes == 0 or handler.get() == nullptr or
        page_handler.get() != nullptr or cache_file != 0)
    {
        Console::error(QObject::tr("Page is already frozen?"));
        return;
    }

    // make sure that the QByteArray takes a copy of the data in d
    // so that we don't get any memory corruption
    auto cached = handler->store(QByteArray(d, nbytes));

    cache_file = cached.first;

    QMutexLocker lkr(&mutex);
    // need to lock before deleting d because fetch() may be called
    delete[] d;
    d = 0;
    page_handler = handler;
    restored_page = cached.second;
}

PageCache PageData::parent() const
{
    return PageCache(c.lock());
}

///////
/////// Implementation of detail::HandleData
///////

/** Construct a HandleData object containing the passed data */
HandleData::HandleData(const QByteArray &data)
    : d(data), offset(0), nbytes(data.size())
{
    if (nbytes == 0)
    {
        state = EMPTY;
    }
    else
    {
        state = DATA_IN_HANDLE;
    }
}

/** Destructor */
HandleData::~HandleData()
{
}

/** Fetch the data from this handle */
QByteArray HandleData::fetch() const
{
    // need to lock because the background thread may be
    // changing the data (moving d to the page)
    QMutexLocker lkr(const_cast<QMutex *>(&mutex));

    if (state == EMPTY or state == DATA_IN_HANDLE)
    {
        return d;
    }

    lkr.unlock();

    // this is a page, so fetch it from here (this call is thread-safe)
    return p.fetch(offset, nbytes);
}

/** Internal function used to set the page containing this data */
void HandleData::setPage(const PageCache::Page &page, int off)
{
    // Need to lock because the background thread will call this,
    // while we may be accessing the data via the thread-safe fetch()
    QMutexLocker lkr(&mutex);

    if (state == EMPTY)
    {
        QString message = QObject::tr("Handle is empty");
        Console::error(message);
        throw SireError::invalid_state(message, CODELOC);
    }
    else if (state == DATA_ON_PAGE)
    {
        QString message = QObject::tr("Handle is already on the page!");
        Console::error(message);
        throw SireError::invalid_state(message, CODELOC);
    }

    p = page;
    offset = off;
    state = DATA_ON_PAGE;
    d = QByteArray();

    lkr.unlock();
}

/** Return the page containing this data. This will be a null
 *  page if the data is not yet on a page
 */
PageCache::Page HandleData::page() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&mutex));
    return p;
}

/** Return the number of bytes of the data held in this handle */
unsigned int HandleData::nBytes() const
{
    return nbytes;
}

/** Return whether or not this handle is valid (actually holds anything) */
bool HandleData::isValid() const
{
    return nbytes > 0;
}

/** Return whether or not this handle is null (does not hold anything) */
bool HandleData::isNull() const
{
    return nbytes == 0;
}

/** Return the cache that this handle is associated with */
PageCache HandleData::parent() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&mutex));
    return p.parent();
}

///////
/////// Implementation of PageCache::Page
///////

/** Construct an empty page */
PageCache::Page::Page()
    : p(nullptr)
{
}

/** Construct a page from the passed data */
PageCache::Page::Page(std::shared_ptr<detail::PageData> data)
    : p(data)
{
}

/** Copy constructor */
PageCache::Page::Page(const PageCache::Page &other)
    : p(other.p)
{
}

/** Destructor */
PageCache::Page::~Page()
{
}

/** Assignment operator */
PageCache::Page &PageCache::Page::operator=(const PageCache::Page &other)
{
    p = other.p;
    return *this;
}

/** Return the type name for this object */
const char *PageCache::Page::typeName()
{
    return "SireBase::PageCache::Page";
}

/** Return the type name for this object */
const char *PageCache::Page::what() const
{
    return PageCache::Page::typeName();
}

/** Return a string representation of this object */
QString PageCache::Page::toString() const
{
    return QString("PageCache::Page(%1 KB used from %2 KB)")
        .arg(this->nBytes() / 1024.0)
        .arg(this->maxBytes() / 1024.0);
}

/** Clone this object */
PageCache::Page *PageCache::Page::clone() const
{
    return new PageCache::Page(*this);
}

/** Assert that this object is valid */
void PageCache::Page::assertValid() const
{
    if (p == nullptr)
    {
        throw SireError::invalid_state(
            QObject::tr("Page object is null"), CODELOC);
    }
}

/** Return whether this page is valid (has some size) */
bool PageCache::Page::isValid() const
{
    return p != nullptr;
}

/** Return whether this page is null (has no size) */
bool PageCache::Page::isNull() const
{
    return p == nullptr;
}

/** Return whether this page is resident in memory */
bool PageCache::Page::isResident() const
{
    if (this->isNull())
    {
        return false;
    }

    return p->isResident();
}

/** Return whether this page is cached to disk */
bool PageCache::Page::isCached() const
{
    if (this->isNull())
    {
        return false;
    }

    return p->isCached();
}

/** Return the number of bytes in this page */
unsigned int PageCache::Page::nBytes() const
{
    if (this->isNull())
    {
        return 0;
    }

    return p->nBytes();
}

/** Return the number of bytes in this page */
unsigned int PageCache::Page::size() const
{
    return this->nBytes();
}

/** Return the maximum number of bytes that can be stored in this page */
unsigned int PageCache::Page::maxBytes() const
{
    if (this->isNull())
    {
        return 0;
    }

    return p->maxBytes();
}

/** Return the parent cache for this page */
PageCache PageCache::Page::parent() const
{
    if (this->isNull())
    {
        return PageCache();
    }

    return p->parent();
}

QByteArray PageCache::Page::fetch(int offset, int n_bytes) const
{
    if (this->isNull())
    {
        throw SireError::invalid_state(
            QObject::tr("Page object is null - you cannot fetch any bytes!"),
            CODELOC);
    }

    return p->fetch(offset, n_bytes);
}

///////
/////// Implementation of PageCache::Handle
///////

/** Construct an empty handle */
PageCache::Handle::Handle()
    : h(nullptr)
{
}

/** Construct a handle from the passed data */
PageCache::Handle::Handle(std::shared_ptr<detail::HandleData> data)
    : h(data)
{
}

/** Copy constructor */
PageCache::Handle::Handle(const PageCache::Handle &other)
    : h(other.h)
{
}

/** Destructor */
PageCache::Handle::~Handle()
{
}

/** Assignment operator */
PageCache::Handle &PageCache::Handle::operator=(const PageCache::Handle &other)
{
    h = other.h;
    return *this;
}

/** Return the type name for this object */
const char *PageCache::Handle::typeName()
{
    return "SireBase::PageCache::Handle";
}

/** Return the type name for this object */
const char *PageCache::Handle::what() const
{
    return PageCache::Handle::typeName();
}

/** Return a string representation of this object */
QString PageCache::Handle::toString() const
{
    const int nbytes = this->nBytes();

    if (nbytes == 0)
    {
        return QString("PageCache::Handle::empty");
    }
    else
    {
        return QString("PageCache::Handle(size = %1 KB)").arg(nbytes / 1024.0);
    }
}

/** Clone this object */
PageCache::Handle *PageCache::Handle::clone() const
{
    return new PageCache::Handle(*this);
}

/** Assert that this object is valid */
void PageCache::Handle::assertValid() const
{
    if (h == nullptr)
    {
        throw SireError::invalid_state(
            QObject::tr("Handle object is null"), CODELOC);
    }
}

/** Return the page on which the data for this handle is placed.
 *  This will be null if the data has not yet been put on a page
 */
PageCache::Page PageCache::Handle::page() const
{
    if (this->isNull())
    {
        return PageCache::Page();
    }

    return Page(h->page());
}

/** Return the data held in this handle */
QByteArray PageCache::Handle::fetch() const
{
    if (this->isNull())
    {
        return QByteArray();
    }

    return h->fetch();
}

/** Return the parent cache for this handle */
PageCache PageCache::Handle::parent() const
{
    if (this->isNull())
    {
        return PageCache();
    }

    return h->parent();
}

/** Return whether this handle is valid (holds data) */
bool PageCache::Handle::isValid() const
{
    return h != nullptr;
}

/** Return whether this handle is null (does not hold data) */
bool PageCache::Handle::isNull() const
{
    return h == nullptr;
}

/** Return the number of bytes in this handle */
unsigned int PageCache::Handle::size() const
{
    return this->nBytes();
}

/** Return the number of bytes in this handle */
unsigned int PageCache::Handle::nBytes() const
{
    if (this->isNull())
    {
        return 0;
    }

    return h->nBytes();
}

/** Clear the data from this handle */
void PageCache::Handle::clear()
{
    h.reset();
}

/** Clear the data in this handle */
void PageCache::Handle::reset()
{
    h.reset();
}

///////
/////// Implementation of PageCache
///////

/** Construct a new page cache with the default
 *  recomended maximum page size */
PageCache::PageCache()
    : d(new CacheData(QString(), CacheData::maxPageSize()))
{
    d->setSelf(d);
}

/** Construct a new page cache with the specified
    recomended maximum page size */
PageCache::PageCache(unsigned int page_size)
    : d(new CacheData(QString(), page_size))
{
    d->setSelf(d);
}

/** Construct a new page cache with specified cache directory
    and recomended maximum page size */
PageCache::PageCache(const QString &cachedir)
    : d(new CacheData(cachedir, CacheData::maxPageSize()))
{
    d->setSelf(d);
}

/** Construct a new page cache with the specified
    recomended maximum page size and cache directory
    template (using QTemporaryDir format) */
PageCache::PageCache(const QString &cachedir, unsigned int page_size)
    : d(new CacheData(cachedir, page_size))
{
    d->setSelf(d);
}

/** Internal constructor used to construct from a CacheData */
PageCache::PageCache(std::shared_ptr<detail::CacheData> data)
    : d(data)
{
}

/** Copy constructor */
PageCache::PageCache(const PageCache &other)
    : d(other.d)
{
}

/** Destructor */
PageCache::~PageCache()
{
}

/** Assignment operator */
PageCache &PageCache::operator=(const PageCache &other)
{
    d = other.d;
    return *this;
}

/** Return the type name for this object */
const char *PageCache::typeName()
{
    return "SireBase::PageCache";
}

/** Return the type name for this object */
const char *PageCache::what() const
{
    return PageCache::typeName();
}

/** Return a string representation of this object */
QString PageCache::toString() const
{
    return QString("PageCache(size = %1 KB)").arg(this->nBytes() / 1024.0);
}

/** Clone this object */
PageCache *PageCache::clone() const
{
    return new PageCache(*this);
}

/** Set the default maximum page cache size for all new created
 *  caches that don't specify it themselves */
void PageCache::setMaxPageSize(unsigned int page_size,
                               bool update_existing)
{
    CacheData::setMaxPageSize(page_size, update_existing);
}

/** Set the maximum number of resident pages per cache */
void PageCache::setMaxResidentPages(unsigned int n_pages)
{
    PageHandler::setMaxResidentPages(n_pages);
}

/** Return the maximum number of resident pages per cache */
unsigned int PageCache::maxResidentPages()
{
    return PageHandler::maxResidentPages();
}

/** Return the current recommend maximum page size */
unsigned int PageCache::maxPageSize()
{
    return CacheData::maxPageSize();
}

/** Assert that this object is valid */
void PageCache::assertValid() const
{
    if (d == nullptr)
    {
        throw SireError::invalid_state(
            QObject::tr("PageCache object is null"), CODELOC);
    }
}

/** Return the statistics for all caches */
QString PageCache::getStatistics()
{
    return CacheData::getStatistics();
}

/** Return the cache directory for this cache */
QString PageCache::cacheDir() const
{
    if (this->isNull())
    {
        return QString();
    }

    return d->cacheDir();
}

/** Return the suggested maximum page size for this cache */
unsigned int PageCache::pageSize() const
{
    if (this->isNull())
    {
        return 0;
    }

    return d->pageSize();
}

/** Return the number of pages in this cache */
unsigned int PageCache::nPages() const
{
    if (this->isNull())
    {
        return 0;
    }

    return d->nPages();
}

/** Return the number of bytes saved in this cache */
unsigned int PageCache::nBytes() const
{
    if (this->isNull())
    {
        return 0;
    }

    return d->nBytes();
}

/** Return the number of bytes saved in this cache */
unsigned int PageCache::size() const
{
    return this->nBytes();
}

/** Return whether or not this cache is valid */
bool PageCache::isValid() const
{
    return d != nullptr;
}

/** Return whether or not this cache is null */
bool PageCache::isNull() const
{
    return d == nullptr;
}

/** Store the data in the cache, returning a handle to the data */
PageCache::Handle PageCache::store(const QByteArray &data)
{
    if (this->isNull())
    {
        // create a cache now
        d = std::make_shared<CacheData>(QString(), 32 * 1024 * 1024);
    }

    return Handle(d->store(data));
}

static QString cache_root_dir = QString();

/** Set the root directory that should be used for all new caches,
 *  when the cache directory is not specified
 */
void PageCache::setRootDirectory(const QString &dir)
{
    QDir d;

    if (not d.mkpath(QDir(dir).absolutePath()))
    {
        throw SireError::io_error(QObject::tr(
                                      "Failed to create cache root directory %1")
                                      .arg(dir),
                                  CODELOC);
    }

    cache_root_dir = QDir(dir).absolutePath();
}

/** Get the root directory that should be used for all new caches,
 *  when the cache directory is not specified
 */
QString PageCache::rootDirectory()
{
    if (cache_root_dir.isEmpty())
    {
        auto env = qEnvironmentVariable("SIRE_PAGECACHE_ROOT", ".");
        PageCache::setRootDirectory(env);
    }

    return cache_root_dir;
}
