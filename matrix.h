#ifndef IOSB_MATRIX_H
#define IOSB_MATRIX_H

#include <cstdlib>
#include <array>
#include <ostream>
#include <complex>
#include <vector>
#include <type_traits>

// - represent a matrix of numerical types (int, long, float, double, complex,....)
// - compute algebraic expressions including + and *, += and *=. Matrix multiplication can be implemented with the most simple algorithm
// - fit in one header file
// - provides a test program which can measure its speed of execution for each examples provided

namespace iosb::storage
{
struct column_major {};
struct row_major    {};
}


namespace iosb::detail {

// \brief expression class for expression templates
//
// \note implements crtp - use of virtual function calls
// 
// \tparam M matrix type for the expression
// \tparam D derived type that can be matrices or generic lambda functions. Must support operator()(std::size_t i)
template<class M, class D>
struct expression
{
	using matrix_type = M;
	decltype(auto) operator()(std::size_t i) const { return static_cast<const D&>(*this)(i); }
	decltype(auto) operator()(std::size_t i)       { return static_cast<      D&>(*this)(i); }
};


// \brief proxy class for encapsulating generic lambdas
// 
// \tparam M matrix type for the lambda function
// \tparam F type of lambda function that is encapsulated
template<class M, class F>
class lambda : public expression <M,lambda<M,F>>
{
public:
	explicit lambda(F const& f)  : expression <M,lambda<M,F>>{}, _f{ f }  {}
	decltype(auto) operator()(std::size_t i) const { return _f(i); }
	decltype(auto) operator()(std::size_t i)       { return _f(i); }
private:
	F _f;
};
// \brief helper function to simply instantiation of lambda proxy class 
template<class M, class F>
auto make_lambda( F&& f ) { return lambda<M,F>(std::forward<F>(f)); }



template<class S, std::size_t M, std::size_t N>
struct index_mapper
{
	static constexpr auto is_column_major = std::is_same<S,::iosb::storage::column_major>::value;
	static constexpr auto is_row_major    = std::is_same<S,::iosb::storage::row_major>::value;
	static_assert(is_column_major || is_row_major, "IndexMapper requires the row or column major tag.");
	static inline constexpr auto wi = is_column_major ? 1 : N;
	static inline constexpr auto wj = is_column_major ? M : 1;
	static inline constexpr auto at(std::size_t i, std::size_t j) { return i * wi + j * wj; }
};


}







namespace iosb 
{

// \brief matrix class
// 
// 
// \tparam E element type
// \tparam M number of rows
// \tparam N number of columns
// \tparam S storage format which can be column_major or row_major
template<class E, std::size_t M, std::size_t N, class S = storage::column_major>
class matrix : public detail::expression<matrix<E,M,N,S>,matrix<E,M,N,S>>
{
public:
	using value_type = E;
	using storage_tag = S;
	using array_type = std::vector<value_type>; // std::array<T,M,N> <- also possible.	
	using reference  = typename array_type::reference;
	using const_reference = typename array_type::const_reference;
	using pointer = typename array_type::pointer;
	using const_pointer = typename array_type::const_pointer;

	
	template<class D>	
	using expression_type  = detail::expression<matrix,D>;
	using base_type = expression_type<matrix>;
	using index_mapper_type = detail::index_mapper<storage_tag,M,N>;
		
	explicit constexpr matrix() : base_type{}, _array(this->size(),value_type{}) {}
	template<class D>
	matrix(expression_type<D> const& other) : base_type{}, _array(this->size(),value_type{}) { this->eval(other);}
	matrix(matrix&& other) : base_type{}, _array{std::move(other._array)} {}	
	
	matrix& operator=(matrix other)
	{
   		swap(*this, other);
		return *this;
	}
	
	friend void swap(matrix& lhs, matrix& rhs)
	{
		using std::swap; 
		swap(lhs._array, rhs._array);
    }	

	template<class D>
	matrix& operator=(expression_type<D> const& other)
	{
		this->eval(other);
		return *this;
	}
	
	~matrix() = default;
	
	inline const_reference operator()(std::size_t i) const { return _array[i]; }
        inline       reference operator()(std::size_t i)       { return _array[i]; }
	
	inline const_reference at(std::size_t ri, std::size_t ci) const { return _array[index_mapper_type::at(ri,ci)]; }
	inline       reference at(std::size_t ri, std::size_t ci)       { return _array[index_mapper_type::at(ri,ci)]; }
	      	
	constexpr auto size() const { return this->rows()*this->cols(); }
	constexpr auto rows() const { return M; }
	constexpr auto cols() const { return N; }
	
	
private:

	template<class D>
	void eval(expression_type<D> const& other)
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
template<class M, class D>
auto make_matrix(detail::expression<M,D> const& other)
{
	return typename detail::expression<M,D>::matrix_type{other};
}
}


// ********* Free Functions ***********

// Matlab Output
template<class E, std::size_t M, std::size_t N, class S>
std::ostream& operator<<(std::ostream& out, iosb::matrix<E,M,N,S> const& m)
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
template<class E, class S, std::size_t M, std::size_t N>
auto operator!(iosb::matrix<E,M,N,S> const& lhs) 
{
	iosb::matrix<E,N,M,S> res{};
	
	for(auto n = 0u; n < res.cols(); ++n)
		for(auto m = 0u; m < res.rows(); ++m)
			res.at(n,m) = lhs.at(m,n);
	
	return res;
}


// Matrix Transpose
template<class E, class S, std::size_t M, std::size_t N, std::size_t Mh, std::size_t Nh>
auto blocked_transpose(iosb::matrix<E,M,N,S> const& lhs) 
{
	
	iosb::matrix<E,N,M,S> res{};
	
	constexpr auto M_ = M/Mh;
	constexpr auto N_ = N/Nh;
	
	if constexpr (std::is_same<S,iosb::storage::column_major>::value)
	{
		
		for(auto m_ = 0u; m_ < M_; ++m_)
			for(auto n_ = 0u; n_ < N_; ++n_)
				for(auto m = m_ * Mh; m < (m_*Mh + Mh); ++m)
					for(auto n = n_ * Nh; n < (n_*Nh + Nh); ++n)
						res.at(n,m) = lhs.at(m,n);						
	}
	else
	{
		
	}
	return res;
}


template<class M,class D>
auto operator!(iosb::detail::expression<M,D> const& lhs) 
{
	return !typename iosb::detail::expression<M,D>::matrix_type{lhs};
}


// Matrix Matrix Multiplikation
template<class E, class S, std::size_t M, std::size_t N, std::size_t K>
auto operator|(iosb::matrix<E,M,K,S> const& lhs, iosb::matrix<E,K,N,S> const& rhs) 
{
	auto res = iosb::matrix<E,M,N,S>{};
	
	if constexpr (std::is_same<S,iosb::storage::column_major>::value)
	{	
		#pragma omp parallel for
		for(auto n = 0u; n < N; ++n)
			for(auto k = 0u; k < K; ++k)
				for(auto m = 0u; m < M; ++m)
					res.at(m,n) += lhs.at(m,k) * rhs.at(k,n);
	}
	else
	{
		#pragma omp parallel for
		for(auto m = 0u; m < M; ++m)
			for(auto k = 0u; k < K; ++k)			
				for(auto n = 0u; n < N; ++n)
					res.at(m,n) += lhs.at(m,k) * rhs.at(k,n);	
	}
	return res;
}


template<class ML, class MR, class L, class R>
auto operator|(iosb::detail::expression<ML,L> const& lhs, iosb::detail::expression<MR,R> const& rhs) 
{
	return typename iosb::detail::expression<ML,L>::matrix_type{lhs}| typename iosb::detail::expression<MR,L>::matrix_type{rhs};
}



// Overloaded arithmetic operators with matrices
template<class M, class L, class R>
auto operator+(iosb::detail::expression<M,L> const& lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M>([&lhs,&rhs](std::size_t i){ return lhs(i) + rhs(i);});
}
template<class M, class L, class R>
auto operator-(iosb::detail::expression<M,L> const& lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M>([&lhs,&rhs](std::size_t i){ return lhs(i) - rhs(i);});
}
template<class M, class L, class R>
auto operator*(iosb::detail::expression<M,L> const& lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M>([&lhs,&rhs](std::size_t i){ return lhs(i) * rhs(i);});
}
template<class M, class L, class R>
auto operator/(iosb::detail::expression<M,L> const& lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M>([&lhs,&rhs](std::size_t i){ return lhs(i) / rhs(i);});
}

// Overloaded Arithmetic Operators with Scalars
template<class M, class R>
auto operator+(typename iosb::detail::expression<M,R>::matrix_type::const_reference lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M>( [&lhs,&rhs](std::size_t i) {return lhs + rhs(i); } );  
}
template<class M, class R>
auto operator-(typename iosb::detail::expression<M,R>::matrix_type::const_reference lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M>( [&lhs,&rhs](std::size_t i) {return lhs - rhs(i); } );
}
template<class M, class R>
auto operator*(typename iosb::detail::expression<M,R>::matrix_type::const_reference lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M> ( [&lhs,&rhs](std::size_t i) {return lhs * rhs(i); } );
}
template<class M, class R>
auto operator/(typename iosb::detail::expression<M,R>::matrix_type::const_reference lhs, iosb::detail::expression<M,R> const& rhs) 
{
	return iosb::detail::make_lambda<M> ( [&lhs,&rhs](std::size_t i) {return lhs / rhs(i); } );
}
template<class M, class L>
auto operator+(iosb::detail::expression<M,L> const& lhs, typename iosb::detail::expression<M,L>::matrix_type::const_reference rhs) 
{
	return iosb::detail::make_lambda<M> ( [&lhs,&rhs](std::size_t i) {return lhs(i) + rhs; } );
}
template<class M, class L>
auto operator-(iosb::detail::expression<M,L> const& lhs, typename iosb::detail::expression<M,L>::matrix_type::const_reference rhs) 
{
	return iosb::detail::make_lambda<M> ( [&lhs,&rhs](std::size_t i) {return lhs(i) - rhs; } );
}
template<class M, class L>
auto operator*(iosb::detail::expression<M,L> const& lhs, typename iosb::detail::expression<M,L>::matrix_type::const_reference rhs) 
{
	return iosb::detail::make_lambda<M> ( [&lhs,&rhs](std::size_t i) {return lhs(i) * rhs; } );
}
template<class M, class L>
auto operator/(iosb::detail::expression<M,L> const& lhs, typename iosb::detail::expression<M,L>::matrix_type::const_reference rhs) 
{
	return iosb::detail::make_lambda<M> ( [&lhs,&rhs](std::size_t i) {return lhs(i) / rhs; } );
}




// Overloaded Assignment Operators
template<class M, class R>
decltype(auto) operator+=(typename iosb::detail::expression<M,R>::matrix_type& lhs, iosb::detail::expression<M,R> const& rhs) { return lhs = lhs + rhs;  }
template<class M, class R>
decltype(auto) operator-=(typename iosb::detail::expression<M,R>::matrix_type& lhs, iosb::detail::expression<M,R> const& rhs) { return lhs = lhs - rhs;  }
template<class M, class R>
decltype(auto) operator*=(typename iosb::detail::expression<M,R>::matrix_type& lhs, iosb::detail::expression<M,R> const& rhs) { return lhs = lhs * rhs;  }
template<class M, class R>
decltype(auto) operator/=(typename iosb::detail::expression<M,R>::matrix_type& lhs, iosb::detail::expression<M,R> const& rhs) { return lhs = lhs / rhs;  }

// Overloaded Assignment Operators
template<class E, class S, std::size_t M, std::size_t N>
decltype(auto) operator+=(iosb::matrix<E,M,N,S>& lhs, typename iosb::matrix<E,M,N,S>::const_reference rhs) { return lhs = lhs + rhs;  }
template<class E, class S, std::size_t M, std::size_t N>
decltype(auto) operator-=(iosb::matrix<E,M,N,S>& lhs, typename iosb::matrix<E,M,N,S>::const_reference rhs) { return lhs = lhs - rhs;  }
template<class E, class S, std::size_t M, std::size_t N>
decltype(auto) operator*=(iosb::matrix<E,M,N,S>& lhs, typename iosb::matrix<E,M,N,S>::const_reference rhs) { return lhs = lhs * rhs;  }
template<class E, class S, std::size_t M, std::size_t N>
decltype(auto) operator/=(iosb::matrix<E,M,N,S>& lhs, typename iosb::matrix<E,M,N,S>::const_reference rhs) { return lhs = lhs / rhs;  }

#endif
