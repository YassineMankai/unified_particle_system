#pragma once

#include "cgp/math/matrix/matrix.hpp"
#include "cgp/containers/buffer_stack/special_types/special_types.hpp"


namespace cgp
{
	using vec9 = buffer_stack<float, 9>;
	using mat39 = matrix_stack<float, 3, 9>;
	using mat9 = matrix_stack<float, 9, 9>;
}
