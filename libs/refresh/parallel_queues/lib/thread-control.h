#ifndef _THREAD_CONTROL
#define _THREAD_CONTROL

#include <mutex>
#include <atomic>
#include <condition_variable>
#include "refresh/parallel_primitives/lib/csemaphore.h"

//#define REFRESH_DISABLE_THREAD_CONTROL

namespace refresh {

	class thread_control
	{
		counting_semaphore_limited semaphore;
	public:
		thread_control(std::ptrdiff_t max_running_at_once) :
			semaphore(max_running_at_once)
		{
			
		}

		std::ptrdiff_t max_running() const
		{
			return semaphore.max_allowed();
		}

		void acquire()
		{
#ifndef REFRESH_DISABLE_THREAD_CONTROL
			semaphore.acquire();
#endif
		}

		void release(std::ptrdiff_t n = 1)
		{
#ifndef REFRESH_DISABLE_THREAD_CONTROL
			semaphore.release(n);
#endif
		}
	};

	//aquires in ctor, releases in dtor
	class thread_control_guard
	{
		thread_control* tc{};
	public:
		thread_control_guard(thread_control& tc): tc(&tc)
		{
			tc.acquire();
		}
		thread_control_guard(const thread_control_guard& rhs) = delete;
		thread_control_guard& operator=(const thread_control_guard& rhs) = delete;

		thread_control_guard(thread_control_guard&& rhs) noexcept : tc(rhs.tc)
		{
			rhs.tc = nullptr;
		}
		thread_control_guard& operator=(thread_control_guard&& rhs) noexcept
		{
			if (this == &rhs) return *this;
			tc = rhs.tc;
			rhs.tc = nullptr;
			return *this;
		}

		~thread_control_guard()
		{
			if (tc)
				tc->release();
		}
	};

	//releases in ctor, acquires in dtor
	class thread_control_guard_reverse
	{
		thread_control* tc{};
	public:
		thread_control_guard_reverse(thread_control& tc) : tc(&tc)
		{
			tc.release();
		}
		thread_control_guard_reverse(const thread_control_guard_reverse& rhs) = delete;
		thread_control_guard_reverse& operator=(const thread_control_guard_reverse& rhs) = delete;

		thread_control_guard_reverse(thread_control_guard_reverse&& rhs) noexcept : tc(rhs.tc)
		{
			rhs.tc = nullptr;
		}
		thread_control_guard_reverse& operator=(thread_control_guard_reverse&& rhs) noexcept
		{
			if (this == &rhs) return *this;
			tc = rhs.tc;
			rhs.tc = nullptr;
			return *this;
		}

		~thread_control_guard_reverse()
		{
			if (tc)
				tc->acquire();
		}
	};

	const uint32_t REFRESH_BUILD_THREAD_CONTROL = 1;
	//very simple class
	//caller side may (should) create treads, but instead of just running them call execute of this class
	//which will just control the number of tasks working at once
	//if it may happen that during execution the thread may be blocked it is worth to call the potentially blocking code
	//in hang_if_needed
	class thread_control_old
	{
		std::atomic<size_t> running{};
		size_t max_running_at_once;
		std::mutex mtx;
		std::condition_variable cv;

		void acquire()
		{
			std::unique_lock<std::mutex> lck(mtx);
			cv.wait(lck, [this] {return running < max_running_at_once; });
			++running;
		}
		void release()
		{
			--this->running;
			cv.notify_one();
		}
	public:
		thread_control_old(size_t max_running_at_once) :
			max_running_at_once(max_running_at_once) {
		}
		//this method may block
		template<typename TASK_T>
		void execute(const TASK_T& task) {
#ifdef REFRESH_DISABLE_THREAD_CONTROL
			task();
#else
			acquire();
			task();
			release();
#endif
		}

		//perform task that may lock current thread (for example when adding to a queue)
		//in this case task should be very cheap because it may be executed beyond the current limit
		//this method may block
		template<typename TASK_T>
		void hang_if_needed(const TASK_T& task)
		{
#ifdef REFRESH_DISABLE_THREAD_CONTROL
			task();
#else
			if (running == max_running_at_once) {
				--running;				//pretend the task is not running
				cv.notify_one(); 		//allow to start some other task (for example consume from queue)
				task(); 				//make this task (for example push to queue which may block) !!! task() must be called before must be before acquire
				acquire(); 				//wait to renew normal work
			}
			else {
				task(); 				//ok there are some free tasks just perform
			}
#endif
		}

		size_t max_running() const {
			return max_running_at_once;
		}
	};
}
#endif // _THREAD_CONTROL
