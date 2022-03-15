#include "custom_shape.hpp"

using namespace cgp;



int shape_structure::getNeighbour(int p_Index, int i,int j){

	int IGrid = (p_Index - relativeLocationsOffset)/width;
	int JGrid= (p_Index - relativeLocationsOffset) % width;
	
	if ((IGrid+i <= height - 1) && (JGrid+j <= width - 1) &&((IGrid+i)>=0) && ((JGrid+j)>=0)) {
		return p_Index + j + i * width;
	}
	
	return -1;

}


int shape_structure::getGlobalIndex(int i, int j) {

	return relativeLocationsOffset + j + i * width;


}



