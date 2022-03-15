#include "helpers_scene.hpp"



int phaseToShapeIndex(int phase)
{
	return std::abs(phase) - 1;
}

int indexShapeToPhase(int indexShape, bool shapeMatching) {
	if (shapeMatching)
		return indexShape + 1;
	else
		return -(indexShape + 1);
}