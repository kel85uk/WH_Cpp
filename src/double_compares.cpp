#include <../include/double_compares.h>

namespace double_compares
{
	bool d_equal_str(realscalar i, realscalar j, realscalar tol){
		return (std::abs((i-j)/i) <= tol) && (std::abs((i-j)/j) <= tol);
	}
	
	bool d_equal_abs_str(realscalar i, realscalar j, realscalar tol){
		return (std::abs((i-j)) <= tol) && (std::abs((i-j)) <= tol);
	}	
	
	bool d_equal_wk(realscalar i, realscalar j, realscalar tol){
		return (std::abs((i-j)/i) <= tol) || (std::abs((i-j)/j) <= tol);
	}	
	
	bool d_equal_str(realscalar i, realscalar j){
		return (std::abs((i-j)/i) <= 100*mach_err) && (std::abs((i-j)/j) <= 100*mach_err);
	}
	
	bool d_equal_abs_str(realscalar i, realscalar j){
		return (std::abs((i-j)) <= 10*mach_err) && (std::abs((i-j)) <= 10*mach_err);
	}	
	
	bool d_equal_wk(realscalar i, realscalar j){
		return (std::abs((i-j)/i) <= 100*mach_err) || (std::abs((i-j)/j) <= 100*mach_err);
	}	
}
