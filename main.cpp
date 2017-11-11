#include "my_matrix.h"

int main() {


	my_matrix A({{1,2,4},{2,1,3},{0,1,1},{-1,3,5}});
	my_matrix B({{1,2,0},{0,1,3},{-1,0,5}});
	A.show(); B.show();
	my_matrix C = A*B;
	C.show();

	my_matrix C1 = B*A;
	C1.show();

	C *= B;
	C.show();
	
	return 0;
}
