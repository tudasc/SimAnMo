#ifndef INCLUDE_xorshf128_TCC_
#define INCLUDE_xorshf128_TCC_

#include <stdint.h>
#include <limits>

namespace blaze_rng
{

template <typename T=uint32_t>
class xorshf128_32bit {
private:
	T x, y, z, w;
	
public:
	typedef T result_type;

	xorshf128_32bit() {
		x=123456789; 
		y=362436069; 
		z=521288629;
		w=88675123;
		for (int i=0; i<10; ++i)
			this->operator ()();
	}
	
	xorshf128_32bit(T _seed) {
		x=123456789; 
		y=362436069; 
		z=521288629;
		w=88675123;
		for (int i=0; i<10; ++i)
			this->operator ()();
	}
	
	result_type operator()() {
	  T t = x ^ (x << 11);
	  x = y; y = z; z = w;
	  w ^= (w >> 19) ^ t ^ (t >> 8);
	  return w;	
	}
	
	static constexpr result_type min () {
		return T(std::numeric_limits<T>::min());
	}
	
	static constexpr result_type max () {
		return T(std::numeric_limits<T>::max());
	}
};

typedef xorshf128_32bit<uint32_t> xorshf128;
} // End namespace


#endif
