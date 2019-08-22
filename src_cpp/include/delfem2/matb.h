#if !defined(MATB_H)
#define MATB_H

#include <assert.h>
#include <vector>

class CTriDiaMat3
{
public:
	// initialize with block size n
	CTriDiaMat3(){ n= 0; v = 0; }
	~CTriDiaMat3(){ delete[] v; }
  void Clear(){ if( n!=0){ delete[] v; } v=0; n=0; }
  void Initialize(int n){
    Clear();
    this->n = n;
    v = new double [n*9*3];
  }
	// clear value
	void SetZero(){ for(unsigned int i=0;i<(n*3-2)*9;i++){ v[i] = 0; } }
	// marge element stiffness matrix to the position (idiv,idiv+1)
	void Marge(unsigned int idiv, double eM[][2][3][3]){
		for(unsigned int i=0;i<36;i++){ v[idiv*27+i] += (&eM[0][0][0][0])[i]; }
	}
	// define fixed boudnary condition
	void FixBC(unsigned int ino, unsigned int idof);
	// execute ILU factorization
	void ILU_Frac();
	// solve matrix
	void Solve(std::vector<double>& res);
private:
	static inline void CalcInvMat3(double a[], double t[] ){
		const double det =
    + a[0]*a[4]*a[8] + a[3]*a[7]*a[2] + a[6]*a[1]*a[5]
    - a[0]*a[7]*a[5] - a[6]*a[4]*a[2] - a[3]*a[1]*a[8];
		const double inv_det = 1.0/det;
		for(int i=0;i<9;i++){ t[i] = a[i]; }
		a[0] = inv_det*(t[4]*t[8]-t[5]*t[7]);
		a[1] = inv_det*(t[2]*t[7]-t[1]*t[8]);
		a[2] = inv_det*(t[1]*t[5]-t[2]*t[4]);
		a[3] = inv_det*(t[5]*t[6]-t[3]*t[8]);
		a[4] = inv_det*(t[0]*t[8]-t[2]*t[6]);
		a[5] = inv_det*(t[2]*t[3]-t[0]*t[5]);
		a[6] = inv_det*(t[3]*t[7]-t[4]*t[6]);
		a[7] = inv_det*(t[1]*t[6]-t[0]*t[7]);
		a[8] = inv_det*(t[0]*t[4]-t[1]*t[3]);
	}
private:
	unsigned int n;
	double* v;
};

#endif
