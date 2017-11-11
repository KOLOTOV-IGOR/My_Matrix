#ifndef MY_MATRIX_H_
#define MY_MATRIX_H_
#include <iostream>
#include <memory>
#include <vector> 
#include <cassert>
#include <omp.h>
using namespace std;

class my_matrix{
private:
	int rows, cols;
	unique_ptr<int[]> M;
public:
	my_matrix() : rows(2), cols(2), M(new int[4]) {
		M[0] = M[2] =  1;
		M[1] = 0; M[3] = -1; 
	}

	my_matrix(const my_matrix &m) {
		rows = m.rows;
		cols = m.cols;
		unique_ptr<int[]> temp(new int[rows*cols]);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				M[i*cols + j] = m.M[i*cols + j];
			}
		}
		M = move(temp);
	}

	my_matrix(int rows, int cols) : rows(rows), cols(cols), M(new int[rows*cols]) {
                for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				M[i*cols + j] = 0;
			}
		}
        }

	my_matrix(int n) : rows(n), cols(n), M(new int[n*n]) {
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				M[i*cols + j] = 0;
			}
		}
	}

	my_matrix(const vector<vector<int>> &a) : rows(a.size()), cols(a[0].size()), M(new int[rows*cols])  {
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				M[i*cols + j] = a[i][j];
			}
		}
	}

	int Rows() const { return rows; }
	int Cols() const { return cols; }

	void identity() {
		if (rows != cols) return;
		for (size_t i = 0; i < rows; i++) {
               	        for (size_t j = 0; j < cols; j++) {
                       	        if (i == j) M[i*cols + j] = 1;
				else M[i*cols + j] = 0;
                       	}
               	}
	}

	void show() {
		for (size_t i = 0; i < rows; i++) {
			cout << "[ ";
                        for (size_t j = 0; j < cols-1; j++) {
                                cout << M[i*cols + j] << ", ";
                        }
			cout << M[i*cols + cols - 1] << " ]" << endl;
                }
		cout << endl;
	}

	my_matrix &operator=(const my_matrix &m) {
		rows = m.rows;
		cols = m.cols;
		unique_ptr<int[]> temp(new int[rows*cols]);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				temp[i*cols + j] = m.M[i*m.cols + j];
			}
		}
		M = move(temp);
		return *this;
	}

	int &operator()(int i, int j) const {
		return M[i*cols + j];
	}

	my_matrix &operator+=(my_matrix &a) {
		if (rows == a.rows && cols == a.cols) {		
			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					M[i*cols + j] += a.M[i*cols + j];
				}
			}
		}
		else { cout << "Error!" << endl;}		
		return *this;
	}

	friend my_matrix &operator+(my_matrix &a, my_matrix &b) {
		a += b;
		return a;
	}

	my_matrix &operator-=(my_matrix &a) {
		if (rows == a.rows && cols == a.cols) {		
			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					M[i*cols + j] -= a.M[i*cols + j];
				}
			}
		}
		else { cout << "Error!" << endl;}		
		return *this;
	}

	friend my_matrix &operator-(my_matrix &a, my_matrix &b) {
		a -= b;
		return a;
	}

	my_matrix &operator*(const int a) {
		for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					M[i*cols + j] *= a;
				}
			}
		return *this;
	}

	friend my_matrix &operator*(const int a, my_matrix &b) {
		for (auto i = 0; i < b.rows; i++) {
				for (auto j = 0; j < b.cols; j++) {
					b(i, j) *= a; 
				} 
		}
		return b;
	}

	friend my_matrix operator*(const my_matrix &a, const my_matrix &b) {
		my_matrix temp(a.rows, b.cols);
		//int n;
		//omp_set_num_threads(11);	
		#pragma omp parallel shared(a, b, temp) //private(n)
		{	
			if (a.cols == b.rows) {
				//n=omp_get_thread_num();
				//printf("%d\n", n);
				#pragma omp for
				for (auto i = 0; i < temp.rows; i++) {
					for (auto j = 0; j < temp.cols; j++) {
						for (auto k = 0; k < a.cols; k++) temp(i, j) += a(i, k)*b(k, j); 
					} 
				}
			}
		}//end of parallel
		return temp;
	}

	friend void operator*=(my_matrix &a, const my_matrix &b) {
		if (a.cols == b.rows) {		
			a = a * b;
		}
		else { cout << "Error!" << endl;}		
	}

	friend my_matrix &Transpose(my_matrix &m) {//for squared matrix
		if (m.rows == m.cols) {
			for (size_t i = 0; i < m.rows; i++) {
               	        	for (size_t j = 0; j < m.cols && i != j; j++) {
			      		swap(m(j, i), m(i, j));			
                       		}
               		}
		}
		else {
			my_matrix temp(m.cols, m.rows);
			for (size_t i = 0; i < m.rows; i++) {
               	        	for (size_t j = 0; j < m.cols; j++) {
			      		temp(j, i) = m(i, j);		
                       		}
               		}
			m = temp;
		}
		return m;
	}

};

#endif




































































