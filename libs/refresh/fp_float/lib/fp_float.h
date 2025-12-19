#pragma once

#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <type_traits>
#include <cmath>
#include <bit>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace refresh {

    // Template class declaration
    template <typename T, std::size_t FRAC_BITS>
    class fp_float {
        static_assert(std::is_integral<T>::value, "T must be an integral type");
        static_assert(std::is_signed<T>::value, "T must be a signed type");
        static_assert(FRAC_BITS > 0 && FRAC_BITS < sizeof(T) * 8, "FRAC_BITS must be positive and less than the number of bits in T");
        static_assert(FRAC_BITS % 2 == 0, "FRAC_BITS must be even");
		static_assert(sizeof(T) == 4 || sizeof(T) == 8, "T must be a 32-bit or 64-bit integer");    
		static_assert(sizeof(T) == 8, "32-bit T is partially supported, needs more testing");    

        constexpr static std::size_t INT_BITS = sizeof(T) * 8;
        constexpr static std::size_t FRAC_BITS_HALF = FRAC_BITS / 2;
        constexpr static T MAX_SHIFTABLE = std::numeric_limits<T>::max() >> FRAC_BITS;
        constexpr static T ONE = T(1) << FRAC_BITS;
        constexpr static T HALF = T(1) << (FRAC_BITS - 1);

        using uT = std::make_unsigned<T>::type;

        T value; // Stored integer value representing the fixed-point number
		static constexpr T error_value = std::numeric_limits<T>::min();

    public:
        bool try_mul(const T& x, const T& y, T& ret) const
        {
#if defined(__GNUC__) || defined(__clang__)
            if (__builtin_mul_overflow(x, y, &ret))
                return false;
#elif defined(_MSC_VER)
            if constexpr (sizeof(T) == 4)
            {
				if (_mul_overflow_i32(x, y, &ret))
					return false;
            }
            else
            {
				if (_mul_overflow_i64(x, y, &ret))
					return false;
            }
#else
            assert(0);
#endif

            return true;
        }

		int no_bits(T x) const
		{
            if(x >= 0)
				return std::bit_width(static_cast<uT>(x));
            else
                return std::bit_width(static_cast<uT>(-x));
		}

        // Friend functions
        template <typename U, std::size_t FB>
        friend fp_float<U, FB> round(const fp_float<U, FB>& x);

        template <typename U, std::size_t FB>
        friend fp_float<U, FB> floor(const fp_float<U, FB>& x);

        template <typename U, std::size_t FB>
        friend fp_float<U, FB> ceil(const fp_float<U, FB>& x);

        template <typename U, std::size_t FB>
        friend fp_float<U, FB> sqrt(const fp_float<U, FB>& x);

        template <typename U, std::size_t FB>
        friend fp_float<U, FB> fast_sqrt(const fp_float<U, FB>& x);

        template <typename U, std::size_t FB>
        friend fp_float<U, FB> abs(const fp_float<U, FB>& x);

        template <typename U, std::size_t FB>
        friend fp_float<U, FB> pow2(const fp_float<U, FB>& x);

        template <typename U, std::size_t FB>
        friend fp_float<U, FB> sub(const fp_float<U, FB>& x, const fp_float<U, FB>& y);

        // Default constructor
        fp_float() : value(0) {}

		fp_float(const fp_float&) = default;
		fp_float(fp_float&&) = default;

        // Constructor from integer types
        template <typename U, typename std::enable_if<std::is_integral<U>::value, int>::type = 0>
        fp_float(U x)
        {
			if (x > std::numeric_limits<T>::max() / ONE)
			{
				value = error_value;
				return;
			}
			else if (x < std::numeric_limits<T>::min() / ONE)
			{
				value = error_value;
				return;
			}

			value = static_cast<T>(x) << FRAC_BITS;
        }

        // Constructor from other fp_float type
        template <std::size_t U_FRAC_BITS>
        fp_float(const fp_float<T, U_FRAC_BITS> &x)
        {
			if (x.is_wrong())
			{
				value = error_value;
				return;
			}

            if (FRAC_BITS == U_FRAC_BITS)
                value = x.raw_value();
            else if (FRAC_BITS > U_FRAC_BITS)
            {
				T max_value = std::numeric_limits<T>::max() >> (FRAC_BITS - U_FRAC_BITS);
				if (x.raw_value() > max_value)
				{
					value = error_value;
					return;
				}
				else if (x.raw_value() < -max_value)
				{
					value = error_value;
					return;
				}   

                value = x.raw_value() << (FRAC_BITS - U_FRAC_BITS);
            }
            else
            {
                value = x.raw_value() >> (U_FRAC_BITS - FRAC_BITS);
            }
        }

        // Constructor from floating-point types
        template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type = 0>
        fp_float(U x)
        {
			if (x > std::numeric_limits<T>::max() / ONE)
			{
				value = error_value;
				return;
			}
			else if (x < std::numeric_limits<T>::min() / ONE)
			{
				value = error_value;
				return;
			}

            value = static_cast<T>(x * ONE);
        }

		fp_float& operator=(const fp_float& other) = default;

		bool is_wrong() const { return value == error_value; }
		
        void make_wrong() { value = error_value; }

		bool is_zero() const { return value == 0; }

        bool is_positive() const { return value > 0; }

        bool is_negative() const { return value < 0; }

        // Operator +=
        fp_float& operator+=(const fp_float& other) {
            if (is_wrong() || other.is_wrong())
            {
				value = error_value;
				return *this;
            }

#if defined(__GNUC__) || defined(__clang__)
            if(__builtin_add_overflow(value, other.value, &value))
				value = error_value;
#elif defined(_MSC_VER)
            if constexpr (sizeof(T) == 4)
            {
				if (_add_overflow_i32(0, value, other.value, &value))
					value = error_value;
            }
            else
            {
				if(_add_overflow_i64(0, value, other.value, &value))
					value = error_value;
            }
#else
            assert(0);
#endif

            return *this;
        }

        // Operator -=
        fp_float& operator-=(const fp_float& other) {
			if (is_wrong() || other.is_wrong())
			{
				value = error_value;
				return *this;
			}

#if defined(__GNUC__) || defined(__clang__)
            if(__builtin_sub_overflow(value, other.value, &value))
				value = error_value;
#elif defined(_MSC_VER)
            if constexpr (sizeof(T) == 4)
            {
				if (_sub_overflow_i32(0, value, other.value, &value))
					value = error_value;
            }
            else
            {
				if(_sub_overflow_i64(0, value, other.value, &value))
					value = error_value;
            }
#else
            assert(0);
#endif

            return *this;
        }

        // Operator +
        fp_float operator+(const fp_float& other) const {
            fp_float result = *this;
            result += other;
            return result;
        }

        // Operator -
        fp_float operator-(const fp_float& other) const {
            fp_float result = *this;
            result -= other;
            return result;
        }

        // Unary operator -
        fp_float operator-() const {
            fp_float result;
            result.value = -value;
            return result;
        }

        // Unary operator +
        fp_float operator+() const {
            return *this;
        }

        // Comparison operators
        bool operator==(const fp_float& other) const {
            return value == other.value;
        }

        bool operator!=(const fp_float& other) const {
            return value != other.value;
        }

        bool operator<(const fp_float& other) const {
            return value < other.value;
        }

        bool operator<=(const fp_float& other) const {
            return value <= other.value;
        }

        bool operator>(const fp_float& other) const {
            return value > other.value;
        }

        bool operator>=(const fp_float& other) const {
            return value >= other.value;
        }

        // Operator *=
        inline fp_float& operator*=(const fp_float& other) {
			if (is_wrong() || other.is_wrong())
			{   
				value = error_value;
				return *this;
			}

			T result;

            if (try_mul(value, other.value, result))
            {
				value = result >> FRAC_BITS;
                return *this;
            }

			int x_bits = no_bits(value);
			int y_bits = no_bits(other.value);
            T x, y;

            if (x_bits < y_bits)
            {
                int x_shift = (int(FRAC_BITS) - (y_bits - x_bits)) / 2;
				if (x_shift < 0)
					x_shift = 0;

				x = value >> x_shift;
				y = other.value >> (int(FRAC_BITS) - x_shift);
            }
            else
            {
				int y_shift = (int(FRAC_BITS) - (x_bits - y_bits)) / 2;
				if (y_shift < 0)
					y_shift = 0;

				x = value >> (int(FRAC_BITS) - y_shift);
				y = other.value >> y_shift;
            }

            if (try_mul(x, y, result))
			{
				value = result;
				return *this;
			}

			value = error_value;

            return *this;
        }

        // Operator *
        inline fp_float operator*(const fp_float& other) const {
            fp_float result = *this;
            result *= other;
            return result;
        }

        // Operator /=
        fp_float& operator/=(const fp_float& other) {
			if (is_wrong() || other.is_wrong() || other.value == 0)
			{
				value = error_value;
				return *this;
			}

            if (value <= MAX_SHIFTABLE && value >= -MAX_SHIFTABLE)
            {
				value <<= FRAC_BITS;
				value /= other.value;
				return *this;
            }

            uT u_numerator;
			uT u_denominator;

            bool negative = false;

			if (value < 0)
			{
				u_numerator = static_cast<uT>(-value);
				negative = !negative;
			}
			else
			{
				u_numerator = value;
			}

			if (other.value < 0)
			{
				u_denominator = static_cast<uT>(-other.value);
				negative = !negative;
			}
			else
			{
				u_denominator = other.value;
			}

			int numerator_shift = INT_BITS - 1 - std::bit_width(u_numerator);

            if (numerator_shift <= (int) FRAC_BITS)
            {
                u_numerator <<= numerator_shift;
                u_denominator >>= FRAC_BITS - numerator_shift;
            }
            else
				u_numerator <<= FRAC_BITS;

			if (u_denominator == 0)
            {
				value = error_value;
				return *this;
			}

			u_numerator /= u_denominator;

			value = negative ? -static_cast<T>(u_numerator) : static_cast<T>(u_numerator);

            return *this;
        }

        // Operator /
        fp_float operator/(const fp_float& other) const {
            fp_float result = *this;
            result /= other;
            return result;
        }

        // Operator >>=
        fp_float& operator>>=(size_t n) {
			value >>= n;
			return *this;
        }

		// Operator >>
		fp_float operator>>(size_t n) const {
			fp_float result = *this;
			result >>= n;
			return result;
		}

        // Operator <<=
        fp_float& operator<<=(size_t n) {
			int x = no_bits(value);
            
            if(x + n < INT_BITS - 1)
    			value <<= n;
            else
				value = error_value;

			return *this;
        }

		// Operator <<
		fp_float operator<<(size_t n) const {
			fp_float result = *this;
			result <<= n;
			return result;
		}

        // Explicit casting operators to integer types
        explicit operator T() const {
            return value >> FRAC_BITS;
        }

        template <typename U, typename std::enable_if<std::is_integral<U>::value, int>::type = 0>
        explicit operator U() const {
            return static_cast<U>(value >> FRAC_BITS);
        }

        // Explicit casting operators to floating-point types
        explicit operator float() const {
            return static_cast<float>(value) / ONE;
        }

        explicit operator double() const {
            return static_cast<double>(value) / ONE;
        }

        // Method to display value (for debugging purposes)
        void display() const {
            std::cout << "fp_float value: " << value << std::endl;
        }

        T raw_value() const { return value; }
    };

    // Function round()
    template <typename T, std::size_t FRAC_BITS>
    fp_float<T, FRAC_BITS> round(const fp_float<T, FRAC_BITS>& x) {
        fp_float<T, FRAC_BITS> result;
        
		if (x.is_wrong())
		{
			result.value = x.error_value;
			return result;
		}

        if (x.value >> (FRAC_BITS - 1) == std::numeric_limits<T>::max() >> (FRAC_BITS - 1))
        {
			result.value = x.error_value;
			return result;
        }

        result.value = x.value + (1 << (FRAC_BITS - 1));
        result.value &= ~((1 << FRAC_BITS) - 1);
        
        return result;
    }

    // Function iround()
    template <typename T, std::size_t FRAC_BITS>
    T int_round(const fp_float<T, FRAC_BITS> &x) 
    {
        assert(!x.is_wrong());

        const T HALF = T(1) << (FRAC_BITS - 1);

        if (x >= 0)
            return ((x.raw_value() >> (FRAC_BITS - 1)) + 1) >> 1;
        else
            return (x.raw_value() + (HALF - 1)) >> FRAC_BITS;
    }

	// Function iround_or_0() - do not crash if x is wrong
    template <typename T, std::size_t FRAC_BITS>
    T int_round_or_0(const fp_float<T, FRAC_BITS> &x) 
    {
        if(x.is_wrong()) [[unlikely]]
			return 0;

        const T HALF = T(1) << (FRAC_BITS - 1);

        if (x >= 0)
            return ((x.raw_value() >> (FRAC_BITS - 1)) + 1) >> 1;
        else
            return (x.raw_value() + (HALF - 1)) >> FRAC_BITS;
    }

    // Function floor()
    template <typename T, std::size_t FRAC_BITS>
    fp_float<T, FRAC_BITS> floor(const fp_float<T, FRAC_BITS>& x) {
        fp_float<T, FRAC_BITS> result;
        
		if (x.is_wrong())
		{
			result.value = x.error_value;
			return result;
		}

        result.value = x.value & ~((1 << FRAC_BITS) - 1);

        return result;
    }

    // Function ceil()
    template <typename T, std::size_t FRAC_BITS>
    fp_float<T, FRAC_BITS> ceil(const fp_float<T, FRAC_BITS>& x) {

        fp_float<T, FRAC_BITS> result;

		if (x.is_wrong())
		{
			result.value = x.error_value;
			return result;
		}

        result.value = x.value;

        if (x.value & ((1 << FRAC_BITS) - 1)) {
            if ((x.value >> FRAC_BITS) == (std::numeric_limits<T>::max() >> FRAC_BITS))
            {
				result.value = x.error_value;
                return result;
            }

            result.value += (1 << FRAC_BITS);
            result.value &= ~((1 << FRAC_BITS) - 1);
        }

        return result;
    }

    // Function sqrt()
    template <typename T, std::size_t FRAC_BITS>
    fp_float<T, FRAC_BITS> sqrt(const fp_float<T, FRAC_BITS>& x) {

//		return fast_sqrt(x);

        fp_float<T, FRAC_BITS> finalResult;
        
		if (x.value < 0)        // no need to check if x is wrong as error_value is negative
        {
			finalResult.value = x.error_value;
			return finalResult;
		}

		using uT = std::make_unsigned<T>::type;

        if(x.value == 0)
		{
			finalResult.value = 0;
			return finalResult;
		}

        uT X = static_cast<uT>(x.value);
		int bw = std::bit_width(X);

        uT s = uT(1) << ((bw + 1) / 2);

        // A few Newton iterations
        // s_{k+1} = floor((s_k + X/s_k)/2)
        for (int i = 0; i < 4; i++) {
            uT d = X / s;
            uT s_next = (s + d) >> 1;
            s = s_next;
        }

        // For sure round up or down
        while (s * s > X)
            s--;
        while ((s + 1) * (s + 1) <= X)
            s++;

        T s_sq = s * s;

        s += (X - s_sq) > (s_sq + 2 * s + 1 - X);

		finalResult.value = static_cast<T>(s << (FRAC_BITS / 2));
		return finalResult;
    }

    // Function sqrt()
    template <typename T, std::size_t FRAC_BITS>
    fp_float<T, FRAC_BITS> fast_sqrt(const fp_float<T, FRAC_BITS>& x) {
        using uT = std::make_unsigned<T>::type;

        fp_float<T, FRAC_BITS> finalResult;
        
		if (x.value <= 0)        // no need to check if x is wrong as error_value is negative
        {
            if(x.value == 0)
				return x;

			finalResult.value = x.error_value;
			return finalResult;
		}

        uT X = static_cast<uT>(x.value);

        uT s = (T) std::sqrt((double)X);

        int bw_thr;
		if constexpr (sizeof(T) == 4)
			bw_thr = 11;
		else
			bw_thr = 25;

        if (s > (uT(1) << bw_thr))
        {
            // A few Newton iterations
            // s_{k+1} = floor((s_k + X/s_k)/2)
//            for (int i = 0; i < 2; i++) {
                uT d = X / s;
                uT s_next = (s + d) >> 1;
                s = s_next;
//            }

            // For sure round up or down
            while (s * s > X)
                s--;
            while ((s + 1) * (s + 1) <= X)
                s++;
        }

        T s_sq = s * s;

        s += (X - s_sq) > (s_sq + 2 * s + 1 - X);

		finalResult.value = static_cast<T>(s << (FRAC_BITS / 2));
		return finalResult;
    }

    // Function abs()
    template <typename T, std::size_t FRAC_BITS>
    fp_float<T, FRAC_BITS> abs(const fp_float<T, FRAC_BITS>& x) {
        fp_float<T, FRAC_BITS> result;

		if (x.is_wrong())
        {
            result.value = x.error_value;
            return result;
        }

        if(x >= 0)
			result.value = x.value;
        else
			result.value = -x.value;

        return result;
    }

    // Function pow2()
    template <typename T, std::size_t FRAC_BITS>
    fp_float<T, FRAC_BITS> pow2(const fp_float<T, FRAC_BITS>& x) {
        fp_float<T, FRAC_BITS> result;

		if (x.is_wrong())
			return x;

        if (x.try_mul(x.value, x.value, result.value))
        {
            result.value >>= FRAC_BITS;
            return result;
        }

		// No need to check if x is wrong as this will be handled by the comparison if shifted max()
        		
        T adjustedValue = x.value >> (FRAC_BITS / 2);

        using uT = std::make_unsigned<T>::type;

        uT y = adjustedValue < 0 ? -adjustedValue : adjustedValue;
        if (y > (std::numeric_limits<T>::max() >> (fp_float<T, FRAC_BITS>::INT_BITS / 2)))
		{
			result.value = x.error_value;
			return result;
		}

        result.value = adjustedValue * adjustedValue;

        return result;
    }

	// Function sub()
    template <typename U, std::size_t FB>
    fp_float<U, FB> sub(const fp_float<U, FB>& x, const fp_float<U, FB>& y)
    {
        fp_float<U, FB> ret;

#if defined(__GNUC__) || defined(__clang__)
        if (__builtin_sub_overflow(x.value, y.value, &ret.value))
            ret.value = ret.error_value;
#elif defined(_MSC_VER)
        if constexpr (sizeof(U) == 4)
        {
			if (_sub_overflow_i32(0, x.value, y.value, &ret.value))
				ret.value = ret.error_value;
        }
        else
        {
			if (_sub_overflow_i64(0, x.value, y.value, &ret.value))
				ret.value = ret.error_value;
        }
#else
        assert(0);
#endif

        return ret;
    }

} // namespace refresh
