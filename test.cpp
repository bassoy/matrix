#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip> 

#include "timer.h"
#include "matrix.h"




int main()
{
	constexpr auto m = 1000u, n = 1000u;
	using matrix = iosb::matrix<double,m,n>;
	using timer  = iosb::timer<iosb::milliseconds>;
	
	auto A = matrix{}, 
	     B = matrix{};
	
	for(auto ri = 0u, j = 0u, k = 0u; ri < A.rows(); ++ri)
		for(auto ci = 0u; ci < A.cols(); ++ci)
			A.at(ri,ci) = k++, B.at(ri,ci) = j++;
	
	timer::tic();
	matrix D = A|B;
	timer::toc();
	std::cout << "Expression[AxB], " << "Time[ms]: "<< timer::elapsed() << "[ms]" << std::endl;
	
	timer::tic();
	matrix E = 2.0 * A - B + B - 4.0*D;
	timer::toc();
	std::cout << "Expression[2.0 * A - B + B- 4.0*D], " << "Time[ms]: "<< timer::elapsed() << "[ms]" << std::endl;
	
	timer::tic();
	matrix F = 2.0 * (A|B) + B + B - 4.0*D + 3.0*E;
	timer::toc();
	std::cout << "Expression[2.0 * AxB + B + B - D- 4.0*D + 3.0*E], " << "Time[ms]: "<< timer::elapsed() << "[ms]" << std::endl;
		
	timer::tic();
	matrix G = (2.0 * A)|(F + B + B -4.0*D + 3.0*E);
	timer::toc();
	std::cout << "Expression[(2.0 * A)x(F + B + B - 4.0*D + 3.0*E)], " << "Time[ms]: "<< timer::elapsed() << "[ms]" << std::endl;
		
	
	if(m<15 && n<15)
	{
		std::ofstream out("check.m");
		out << std::setprecision(15);
		out << "A=" << A << std::endl;
		out << "B=" << B << std::endl;
		out << "D=" << D << std::endl;
		out << "E=" << E << std::endl;
		out << "F=" << F << std::endl;
		out << "G=" << G << std::endl;
		out << "Dref = A*B;" << std::endl;
		out << "Eref = 2.*A.-B.+B.-4.*Dref;" << std::endl;
		out << "Fref = 2.*A*B.+B.+B.-4.*Dref.+3.*Eref;" << std::endl;
		out << "Gref = (2.*A)*(Fref.+B.+B.-4.*Dref.+3.*Eref);" << std::endl << std::endl;
		out << "printf('MaxAbsErr(D) = %f\\n',max(D(:)-Dref(:)));" << std::endl;
		out << "printf('MaxAbsErr(E) = %f\\n',max(E(:)-Eref(:)));" << std::endl;
		out << "printf('MaxAbsErr(F) = %f\\n',max(F(:)-Fref(:)));" << std::endl;
		out << "printf('MaxAbsErr(G) = %f\\n',max(G(:)-Gref(:)));" << std::endl;
	}
}
