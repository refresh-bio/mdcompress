#pragma once

#include <atomic>
#include <cstddef>

namespace refresh
{
    class counting_semaphore_limited 
    {
        std::ptrdiff_t init_counter;
        std::atomic<std::ptrdiff_t> counter;
        std::atomic<std::ptrdiff_t> no_waiting;

        void local_wait() noexcept
        {
            no_waiting.fetch_add(1);

            auto current = counter.load();
            if (!current)
                atomic_wait(&counter, current);

            no_waiting.fetch_sub(1, std::memory_order_relaxed);
        }

    public:
        constexpr explicit counting_semaphore_limited(const ptrdiff_t desired) noexcept
            : init_counter(desired), counter(desired)
        {
        }

        counting_semaphore_limited(const counting_semaphore_limited&) = delete;
        counting_semaphore_limited& operator=(const counting_semaphore_limited&) = delete;

		std::ptrdiff_t max_allowed() const noexcept
		{
			return init_counter;
		}

        void release(std::ptrdiff_t update = 1) noexcept
        {
            if (update == 0)
                return;

            counter++;
            const auto waiting_upper_bound = no_waiting.load();

            if (waiting_upper_bound == 0)
                ;
            else if (waiting_upper_bound <= update)
                counter.notify_all();
            else
                while(update--)
                    counter.notify_one();
        }

        void acquire() noexcept
        {
            auto current = counter.load(std::memory_order_relaxed);
            
            while(true)
            {
                while (current == 0) 
                {
                    local_wait();
                    current = counter.load(std::memory_order_relaxed);
                }

                if (counter.compare_exchange_weak(current, current - 1))
                    return;
            }
        }

        bool try_acquire() noexcept
        {
            auto current = counter.load();
            
            if (current == 0)
                return false;

            return counter.compare_exchange_weak(current, current - 1);
        }
    };
}