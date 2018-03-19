#ifndef IOSB_MATRIX_H
#define IOSB_MATRIX_H

#include <cstdlib>
#include <array>
#include <ostream>
#include <complex>
#include <vector>

// - represent a matrix of numerical types (int, long, float, double, complex,....)
// - compute algebraic expressions including + and *, += and *=. Matrix multiplication can be implemented with the most simple algorithm
// - fit in one header file
// - provides a test program which can measure its speed of execution for each examples provided

namespace iosb {
namespace detail {

// \brief expression class for expression templates
//
// \note implements crtp - use of virtual function calls
// 
// \tparam T element type of matrices and scalars of the expression
// \tparam D derived type that can be matrices or generic lambda functions. Must support operator()(std::size_t i)
template<class T, std::size_t M, std::size_t N, class D>
struct expression
{
	decltype(auto) operator()(std::size_t i) const { return static_cast<const D&>(*this)(i); }
	decltype(auto) operator()(std::size_t i)       { return static_cast<      D&>(*this)(i); }
};
// \brief proxy class for encapsulating generic lambdas
// 
// \tparam T element type of matrices and scalars of the expression
// \tparam F type of lambda function that is encapsulated
template<class T, std::size_t M, std::size_t N, class F>
class lambda : public expression <T,M,N,lambda<T,M,N,F> >
{
public:
	using value_type = T;
	using expr_type  = expression <T,M,N,lambda<T,M,N,F> >;
	explicit lambda(F const& f)  : expr_type{}, _f{ f }  {}
	decltype(auto) operator()(std::size_t i) const { return _f(i); }
	decltype(auto) operator()(std::size_t i)       { return _f(i); }
private:
	F _f;
};
// \brief helper function to simply instantiation of lambda proxy class 
template<class T, std::size_t M, std::size_t N,  class F>
auto make_lambda( F&& f ) { return lambda<T,M,N,F>(std::forward<F>(f)); }
}
}


namespace iosb 
{
// \brief matrix class
// 
// \note here we have column-major storage format
// 
// \tparam T element type
// \tparam M number of rows
// \tparam N number of columns
template<class T, std::size_t M, std::size_t N>
class matrix : public detail::expression<T,M,N,matrix<T,M,N> >
{
public:
	using array_type = std::vector<T>; // std::array<T,M,N> <- also possible.
	using expr_type  = detail::expression<T,M,N,matrix<T,M,N>>;
	using value_type = typename array_type::value_type;
	using reference  = typename array_type::reference;
	using const_reference = typename array_type::const_reference;
	using pointer = typename array_type::pointer;
	using const_pointer = typename array_type::const_pointer;
		
	explicit constexpr matrix() : expr_type{}, _array(M*N,T{}) {}
	template<class D>
	matrix(detail::expression<T,M,N,D> const& other) : expr_type{}, _array(M*N,T{}) { this->eval(other);}
	matrix(matrix&& other) : expr_type{}, _array{std::move(other._array)} {}	
	
	matrix& operator=(matrix other)
	{
   		other.swap (*this);
		return *this;
	}
	
	friend void swap(matrix& lhs, matrix& rhs){
		using std::swap; 
		swap(lhs._array, rhs._array);
    }	

	template<class D>
	matrix& operator=(detail::expression<T,M,N,D> const& other)
	{
		this->eval(other);
		return *this;
	}
	
	~matrix() = default;
	
	const_reference operator()(std::size_t i) const { return _array[i]; }
          reference operator()(std::size_t i)       { return _array[i]; }
	
	const_reference at(std::size_t ri, std::size_t ci) const { return _array[ri + this->rows() * ci]; }
	      reference at(std::size_t ri, std::size_t ci)       { return _array[ri + this->rows() * ci]; }
	      	
	constexpr auto size() const { return M*N; }
	constexpr auto rows() const { return M; }
	constexpr auto cols() const { return N; }
	
	
private:

	template<class D>
	void eval(detail::expression<T,M,N,D> const& other)
	{
		#pragma omp parallel for
		for(auto i = 0u; i < this->size(); ++i)
			_array[i] = other(i);
	}
	
	array_type _array;
};
}


namespace iosb 
{
template<class T, std::size_t M, std::size_t N, class D>
auto make_matrix(detail::expression<T,M,N,D> const& other)
{
	return matrix<T,M,N>{other};
}
}


// ********* Free Functions ***********

// Matlab Output
template<class T, std::size_t M, std::size_t N>
std::ostream& operator<<(std::ostream& out, iosb::matrix<T,M,N> const& m)
{
	out << "[ ... " << std::endl;
	for(auto ri = 0u; ri < m.rows(); ++ri){
		for(auto ci = 0u; ci < m.cols(); ++ci)
			out << m.at(ri, ci) << " ";
		out << "; ... " << std::endl;
	}
	out << "];" << std::endl;
	return out;
}

// Matrix Transpose
template<class T, std::size_t M, std::size_t N>
auto operator!(iosb::matrix<T,M,N> const& lhs) 
{
	iosb::matrix<T,N,M> res{};
	
	for(auto n = 0u; n < N; ++n)
		for(auto m = 0u; m < M; ++m)
			res.at(n,m) = lhs.at(m,n);
	
	return res;
}

// Matrix Matrix Multiplikation
template<class T, std::size_t M, std::size_t N, std::size_t K>
auto operator|(iosb::matrix<T,M,K> const& lhs, iosb::matrix<T,K,N> const& rhs) 
{
	iosb::matrix<T,M,N> res{};
	
	#pragma omp parallel for
	for(auto n = 0u; n < N; ++n)
		for(auto k = 0u; k < K; ++k)
			for(auto m = 0u; m < M; ++m)
				res.at(m,n) += lhs.at(m,k) * rhs.at(k,n);
	return res;
}

template<class T, std::size_t M, std::size_t N, std::size_t K, class L, class R>
auto operator|(iosb::detail::expression<T,M,K,L> const& lhs, iosb::detail::expression<T,K,N,R> const& rhs) 
{
	return iosb::matrix<T,M,K>{lhs}|iosb::matrix<T,K,N>{rhs};
}



// Overloaded arithmetic operators with matrices
template<class T, std::size_t M, std::size_t N, class L, class R>
auto operator+(iosb::detail::expression<T,M,N,L> const& lhs,iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N>([&lhs,&rhs](std::size_t i){ return lhs(i) + rhs(i);});
}

template<class T, std::size_t M, std::size_t N, class L, class R>
auto operator-(iosb::detail::expression<T,M,N,L> const& lhs,iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N>([&lhs,&rhs](std::size_t i){ return lhs(i) - rhs(i);});
}

template<class T, std::size_t M, std::size_t N, class L, class R>
auto operator*(iosb::detail::expression<T,M,N,L> const& lhs,iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N>([&lhs,&rhs](std::size_t i){ return lhs(i) * rhs(i);});
}

template<class T, std::size_t M, std::size_t N, class L, class R>
auto operator/(iosb::detail::expression<T,M,N,L> const& lhs,iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N>([&lhs,&rhs](std::size_t i){ return lhs(i) / rhs(i);});
}

// Overloaded Arithmetic Operators with Scalars
template<class T, std::size_t M, std::size_t N, class R>
auto operator+(T const& lhs, iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs + rhs(i); } );  
}

template<class T, std::size_t M, std::size_t N, class R>
auto operator-(T const& lhs, iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs - rhs(i); } );
}

template<class T, std::size_t M, std::size_t N, class R>
auto operator*(T const& lhs, iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs * rhs(i); } );
}

template<class T, std::size_t M, std::size_t N, class R>
auto operator/(T const& lhs, iosb::detail::expression<T,M,N,R> const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs / rhs(i); } );
}

template<class T, std::size_t M, std::size_t N, class L>
auto operator+(iosb::detail::expression<T,M,N,L> const& lhs, T const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs(i) + rhs; } );
}

template<class T, std::size_t M, std::size_t N, class L>
auto operator-(iosb::detail::expression<T,M,N,L> const& lhs, T const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs(i) - rhs; } );
}

template<class T, std::size_t M, std::size_t N, class L>
auto operator*(iosb::detail::expression<T,M,N,L> const& lhs, T const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs(i) * rhs; } );
}

template<class T, std::size_t M, std::size_t N, class L>
auto operator/(iosb::detail::expression<T,M,N,L> const& lhs, T const& rhs) 
{
	return iosb::detail::make_lambda<T,M,N> ( [&lhs,&rhs](std::size_t i) {return lhs(i) / rhs; } );
}


// Overloaded Assignment Operators
template<class T, std::size_t M, std::size_t N, class R>
decltype(auto) operator+=(iosb::matrix<T,M,N>& lhs, iosb::detail::expression<T,M,N,R> const& rhs) { return lhs = lhs + rhs;  }
template<class T, std::size_t M, std::size_t N, class R>
decltype(auto) operator-=(iosb::matrix<T,M,N>& lhs, iosb::detail::expression<T,M,N,R> const& rhs) { return lhs = lhs - rhs;  }
template<class T, std::size_t M, std::size_t N, class R>
decltype(auto) operator*=(iosb::matrix<T,M,N>& lhs, iosb::detail::expression<T,M,N,R> const& rhs) { return lhs = lhs * rhs;  }
template<class T, std::size_t M, std::size_t N, class R>
decltype(auto) operator/=(iosb::matrix<T,M,N>& lhs, iosb::detail::expression<T,M,N,R> const& rhs) { return lhs = lhs / rhs;  }

// Overloaded Assignment Operators
template<class T, std::size_t M, std::size_t N>
decltype(auto) operator+=(iosb::matrix<T,M,N>& lhs, T const& rhs) { return lhs = lhs + rhs;  }
template<class T, std::size_t M, std::size_t N, class R>
decltype(auto) operator-=(iosb::matrix<T,M,N>& lhs, T const& rhs) { return lhs = lhs - rhs;  }
template<class T, std::size_t M, std::size_t N, class R>
decltype(auto) operator*=(iosb::matrix<T,M,N>& lhs, T const& rhs) { return lhs = lhs * rhs;  }
template<class T, std::size_t M, std::size_t N, class R>
decltype(auto) operator/=(iosb::matrix<T,M,N>& lhs, T const& rhs) { return lhs = lhs / rhs;  }

#endif
