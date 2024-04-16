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

#include <QDir>
#include <QThread>
#include <QQueue>
#include <QMutex>

#include <boost/noncopyable.hpp>

namespace SireBase
{
    namespace detail
    {
        class CacheData : public QThread
        {
        public:
            CacheData(QString cachedir, int page_size);
            ~CacheData();

            static QString getStatistics();

            QString cacheDir() const;

            PageCache::Handle cache(const QByteArray &data);

            int pageSize() const;

            int nPages() const;
            int nBytes() const;

            void registerCache(const std::shared_ptr<CacheData> &cache);

        protected:
            void run();

        private:
            void enqueue(const std::shared_ptr<HandleData> &handle);

            static QMutex caches_mutex;
            static QList<std::weak_ptr<CacheData>> caches;

            QMutex queue_mutex;
            QQueue<std::shared_ptr<HandleData>> queue;

            QMutex page_mutex;
            std::weak_ptr<PageData> current_page;
            QList<std::weak_ptr<PageData>> pages;

            std::weak_ptr<CacheData> self;

            QString cache_dir;
            int page_size;

            bool exiting;
        };

        class PageData : public boost::noncopyable
        {
        public:
            PageData(int max_size, const std::shared_ptr<CacheData> &cache);
            ~PageData();

            int maxBytes() const;
            int nBytes() const;
            int bytesRemaining() const;

            bool isResident() const;
            bool isCached() const;

            int cache(const QByteArray &data);

            QByteArray fetch(int offset, int n_bytes) const;

            PageCache parent() const;

        private:
            std::shared_ptr<CacheData> c;

            int max_bytes;
            int nbytes;

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

            int nBytes() const;

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

QMutex CacheData::caches_mutex;
QList<std::weak_ptr<CacheData>> CacheData::caches;

void CacheData::registerCache(const std::shared_ptr<CacheData> &cache)
{
    QMutexLocker lkr(&caches_mutex);
    caches.append(cache);
    self = cache;
}

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

CacheData::CacheData(QString c, int p)
    : cache_dir(c), page_size(p), exiting(false)
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
    this->requestInterruption();

    if (QThread::currentThread() != this)
    {
        this->wait();
    }
    else
    {
        qWarning() << "CacheData is deleting itself!" << this->cacheDir();
    }
}

PageCache::Handle CacheData::cache(const QByteArray &data)
{
    auto handle = std::make_shared<HandleData>(data);
    this->enqueue(handle);
    return PageCache::Handle(handle);
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
    QMutexLocker lkr(const_cast<QMutex *>(&page_mutex));
    int npages = 0;

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

int CacheData::nBytes() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&page_mutex));

    int nbytes = 0;

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

    if (not this->isRunning())
    {
        this->start();
    }
}

void CacheData::run()
{
    // get hold of a pointer to self, so that we aren't
    // deleted while we are running
    auto locked_self = self.lock();

    int empty_count = 0;

    while (true)
    {
        if (this->isInterruptionRequested())
        {
            // stop what we are doing
            return;
        }

        std::shared_ptr<HandleData> handle;
        bool have_item = false;

        // pull an item off the queue
        {
            QMutexLocker lkr(&queue_mutex);

            if (this->isInterruptionRequested())
            {
                // stop what we are doing
                return;
            }

            while (not queue.isEmpty())
            {
                handle = queue.dequeue();

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

            const int n_bytes = data.size();

            if (n_bytes >= page_size)
            {
                // this is bigger than a page, so needs to have its
                // own page!
                auto page = std::make_shared<PageData>(n_bytes, this->self.lock());
                QMutexLocker lkr(&page_mutex);
                this->pages.append(page);
                lkr.unlock();
                auto offset = page->cache(data);
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
                    auto offset = current->cache(data);
                    handle->setPage(PageCache::Page(current), offset);
                }
                else
                {
                    // add the current page to the list of old pages
                    QMutexLocker lkr(&page_mutex);
                    this->pages.append(current_page);

                    // we need to create a new page
                    current = std::make_shared<PageData>(page_size, this->self.lock());
                    current_page = current;
                    lkr.unlock();

                    auto offset = current->cache(data);
                    handle->setPage(PageCache::Page(current), offset);
                }
            }

            empty_count = 0;
        }

        if (this->isInterruptionRequested())
        {
            // stop what we are doing
            return;
        }

        bool is_empty = true;

        if (have_item)
        {
            // check to see if there is anything else to process
            QMutexLocker lkr(&queue_mutex);
            is_empty = queue.isEmpty();
            empty_count += 1;

            if (empty_count > 10)
            {
                // we have been idle for a while, so we can stop
                this->exiting = true;
                return;
            }

            lkr.unlock();
        }

        if (is_empty)
        {
            if (this->isInterruptionRequested())
            {
                // stop what we are doing
                return;
            }

            // sleep for a bit
            this->msleep(100);
        }
    }
}

///////
/////// Implementation of detail::PageData
///////

PageData::PageData(int max_size, const std::shared_ptr<CacheData> &cache)
    : c(cache), max_bytes(max_size), nbytes(0), d(0)
{
    if (max_bytes < 1024)
    {
        throw SireError::invalid_arg(
            QObject::tr("Page size must be greater than 1024!"), CODELOC);
    }
    else if (max_bytes > 128 * 1024 * 1024)
    {
        throw SireError::invalid_arg(
            QObject::tr("Page size must be less than 128 MB!"), CODELOC);
    }

    d = new char[max_bytes];
}

PageData::~PageData()
{
    delete[] d;
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

int PageData::cache(const QByteArray &data)
{
    if (data.size() > this->bytesRemaining())
    {
        throw SireError::invalid_arg(
            QObject::tr("Data is too large to fit on this page!"), CODELOC);
    }

    std::memcpy(d + nbytes, data.constData(), data.size());

    int offset = nbytes;
    nbytes += data.size();

    return offset;
}

QByteArray PageData::fetch(int offset, int n_bytes) const
{
    if (offset + n_bytes > nbytes)
    {
        throw SireError::invalid_arg(
            QObject::tr("Data is too large to fit on this page!"), CODELOC);
    }

    return QByteArray(d + offset, n_bytes);
}

PageCache PageData::parent() const
{
    return PageCache(c);
}

///////
/////// Implementation of detail::HandleData
///////

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

HandleData::~HandleData()
{
}

QByteArray HandleData::fetch() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&mutex));

    if (state == EMPTY or state == DATA_IN_HANDLE)
    {
        return d;
    }

    lkr.unlock();

    return p.fetch(offset, nbytes);
}

void HandleData::setPage(const PageCache::Page &page, int off)
{
    QMutexLocker lkr(&mutex);

    if (state == EMPTY)
    {
        throw SireError::invalid_state(
            QObject::tr("Handle is empty"), CODELOC);
    }
    else if (state == DATA_ON_PAGE)
    {
        throw SireError::invalid_state(
            QObject::tr("Handle is already on the page!"), CODELOC);
    }

    p = page;
    offset = off;
    state = DATA_ON_PAGE;
    d = QByteArray();

    lkr.unlock();
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
    return QString("PageCache::Page(%1 KB used from %2 KB)")
        .arg(this->nBytes() / 1024.0)
        .arg(this->maxBytes() / 1024.0);
}

PageCache::Page *PageCache::Page::clone() const
{
    return new PageCache::Page(*this);
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

int PageCache::Page::maxBytes() const
{
    assertValid();
    return p->maxBytes();
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

PageCache::Handle *PageCache::Handle::clone() const
{
    return new PageCache::Handle(*this);
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
    d->registerCache(d);
}

PageCache::PageCache(const QString &cachedir, int page_size)
    : d(new CacheData(cachedir, page_size))
{
    d->registerCache(d);
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
    return QString("PageCache(size = %1 KB)").arg(this->nBytes() / 1024.0);
}

PageCache *PageCache::clone() const
{
    return new PageCache(*this);
}

void PageCache::assertValid() const
{
    if (d == nullptr)
    {
        throw SireError::invalid_state(
            QObject::tr("PageCache object is null"), CODELOC);
    }
}

QString PageCache::getStatistics()
{
    return CacheData::getStatistics();
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

PageCache::Handle PageCache::store(const QByteArray &data)
{
    return this->cache(data);
}
