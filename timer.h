#ifndef IOSB_TIMER_H
#define IOSB_TIMER_H


namespace iosb
{

using nanoseconds  = std::chrono::duration<double, std::nano>;
using microseconds = std::chrono::duration<double, std::micro>;
using milliseconds = std::chrono::duration<double, std::milli>;
using seconds      = std::chrono::duration<double, std::ratio<1,1>>;
template<class D>
class timer
{
public:
	using clock = std::chrono::high_resolution_clock;
	using point = typename clock::time_point;
	
	static void tic() { _start = clock::now(); }
	static void toc() { _stop  = clock::now(); }

	static auto elapsed() { return std::chrono::duration_cast<D>( _stop - _start ); }

	template<class O>
	static auto elapsed() { return std::chrono::duration_cast<O>( _stop - _start ); }

	static inline point start() { return _start; }
	static inline point stop() { return _stop; }
private:
	static point _start;
	static point _stop;
};
template<class D> typename timer<D>::point timer<D>::_start;
template<class D> typename timer<D>::point timer<D>::_stop;
}
template<class T>
std::ostream& operator<< (std::ostream & out, std::chrono::duration<double,T> const& d) 
{ out << d.count(); return out; }


#endif
