#pragma once

// CPU-only BoxData storage reuse. Accelerator builds retain their existing
// allocation/free semantics, including any synchronization in those paths.
#ifndef PROTO_ACCEL
#include <cstdlib>
#include <memory>
#include <mutex>
#include <vector>

#ifndef PR_HOST_BUFFER_CACHE_BYTES
#define PR_HOST_BUFFER_CACHE_BYTES (64ULL * 1024 * 1024)
#endif

namespace Proto { namespace detail {

class HostBufferCache
{
public:
    static constexpr size_t minimumBytes = 64 * 1024;

    explicit HostBufferCache(size_t a_limit = PR_HOST_BUFFER_CACHE_BYTES)
      : m_limit(a_limit), m_cachedBytes(0) {}

    ~HostBufferCache()
    {
        for (const auto& entry : m_buffers) { std::free(entry.pointer); }
    }

    HostBufferCache(const HostBufferCache&) = delete;
    HostBufferCache& operator=(const HostBufferCache&) = delete;

    void* allocate(size_t a_bytes)
    {
        if (eligible(a_bytes))
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            // Recent exact-size matches tend to be the next patch's temporaries.
            for (size_t ii = m_buffers.size(); ii > 0; --ii)
            {
                if (m_buffers[ii - 1].bytes == a_bytes)
                {
                    void* pointer = m_buffers[ii - 1].pointer;
                    m_buffers.erase(m_buffers.begin() + ii - 1);
                    m_cachedBytes -= a_bytes;
                    return pointer;
                }
            }
        }
        return std::malloc(a_bytes);
    }

    void release(void* a_pointer, size_t a_bytes) noexcept
    {
        if (!a_pointer) { return; }
        if (!eligible(a_bytes)) { std::free(a_pointer); return; }
        try
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            // Evict the oldest idle buffers; the limit bounds retained bytes,
            // independently of how much live application storage exists.
            while (m_cachedBytes > m_limit - a_bytes)
            {
                auto entry = m_buffers.front();
                m_buffers.erase(m_buffers.begin());
                m_cachedBytes -= entry.bytes;
                std::free(entry.pointer);
            }
            m_buffers.push_back({a_pointer, a_bytes});
            m_cachedBytes += a_bytes;
        }
        catch (...)
        {
            // Failure to grow cache bookkeeping must not throw from a deleter.
            std::free(a_pointer);
        }
    }

    size_t cachedBytes() const
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_cachedBytes;
    }

private:
    bool eligible(size_t a_bytes) const
    {
        return a_bytes >= minimumBytes && a_bytes <= m_limit;
    }

    struct Entry { void* pointer; size_t bytes; };
    const size_t m_limit;
    size_t m_cachedBytes;
    std::vector<Entry> m_buffers;
    mutable std::mutex m_mutex;
};

// Keep a tiny control shell alive through static destruction. The exit callback
// releases its cache ownership; late static allocations then use malloc/free.
// Existing buffer owners can safely outlive the callback through shared ownership.
struct HostBufferCacheState
{
    std::mutex mutex;
    std::shared_ptr<HostBufferCache> cache = std::make_shared<HostBufferCache>();
};

inline void shutdownHostBufferCache();

inline HostBufferCacheState& hostBufferCacheState()
{
    static auto* state = []() {
        auto* control = new HostBufferCacheState;
        if (std::atexit(shutdownHostBufferCache) != 0) { control->cache.reset(); }
        return control;
    }();
    return *state;
}

inline void shutdownHostBufferCache()
{
    auto& state = hostBufferCacheState();
    std::lock_guard<std::mutex> lock(state.mutex);
    state.cache.reset();
}

inline std::shared_ptr<HostBufferCache> hostBufferCache()
{
    auto& state = hostBufferCacheState();
    std::lock_guard<std::mutex> lock(state.mutex);
    return state.cache;
}

template<typename T>
inline std::shared_ptr<T> allocateHostBoxData(size_t a_bytes)
{
    if (a_bytes < HostBufferCache::minimumBytes ||
        a_bytes > PR_HOST_BUFFER_CACHE_BYTES)
    {
        return std::shared_ptr<T>(static_cast<T*>(std::malloc(a_bytes)),
            [](T* pointer) { std::free(pointer); });
    }
    auto cache = hostBufferCache();
    if (!cache)
    {
        return std::shared_ptr<T>(static_cast<T*>(std::malloc(a_bytes)),
            [](T* pointer) { std::free(pointer); });
    }
    auto pointer = static_cast<T*>(cache->allocate(a_bytes));
    // Holding the cache here also makes destruction of static BoxData owners
    // safe: idle storage is freed only after the last owner releases the cache.
    return std::shared_ptr<T>(pointer, [cache, a_bytes](T* data)
        { cache->release(data, a_bytes); });
}

}} // namespace Proto::detail
#endif
