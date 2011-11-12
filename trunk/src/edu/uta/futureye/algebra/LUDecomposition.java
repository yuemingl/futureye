package edu.uta.futureye.algebra;

import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;

public class LUDecomposition {

	/*
	 * A=L*U 
	 * L is a lower triangular matrix with ones on its diagonal
	 * U is an upper triangular matrix
	 *  
	 * Note: it is without permutation of L (this function is just for test only)
	 */
	public static void LU(Matrix A, Matrix L, Matrix U) {
		int N = A.getRowDim();
		
		for(int n=1; n<=N; n++) {
			for(int j=n; j<=N; j++) {
				double sum1 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum1 += L.get(n, k)*U.get(k, j);
				}
				double anj = A.get(n, j);
				U.set(n, j, anj-sum1);
			}
			for(int i=n+1; i<=N; i++) {
				double sum2 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum2 += L.get(i, k)*U.get(k, n);
				}
				double ain = A.get(i, n);
				L.set(i, n, (ain-sum2)/U.get(n, n));
			}
			L.set(n, n, 1.0);
		}
	}
	
	public static void LU(FullMatrix A, FullMatrix L, FullMatrix U) {
		int N = A.getRowDim();
		double[][] dA = A.data;
		double[][] dL = L.data;
		double[][] dU = U.data;
		
		for(int n=0; n<=N-1; n++) {
			for(int j=n; j<N; j++) {
				double sum1 = 0.0;
				for(int k=0; k<=n-1; k++) {
					sum1 += dL[n][k]*dU[k][j];
				}
				double anj = dA[n][j];
				dU[n][j] = anj-sum1;
			}
			for(int i=n+1; i<N; i++) {
				double sum2 = 0.0;
				for(int k=0; k<=n-1; k++) {
					sum2 += dL[i][k]*dU[k][n];
				}
				double ain = dA[i][n];
				dL[i][n] = (ain-sum2)/dU[n][n];
			}
			dL[n][n] = 1.0;
		}
	}
	
	/*
	 * A=P*L*U
	 * 
	 * P is a permutation matrix
	 * L is a lower triangular matrix with ones on its diagonal
	 * U is an upper triangular matrix
	 * 
	 */
	public static void LU(Matrix AA, Matrix L, Matrix U, Matrix P) {
		Matrix A = AA.copy();
		int N = A.getRowDim();
		int[] VP = new int[N+1];
		for(int n=1; n<=N; n++) {
			VP[n] = n;
		}
		for(int n=1; n<=N; n++) {
			//判断Unn的值是否为0，为0的话进行行对换，直到找到Unn不为0的行
			for(int nn=n; nn<=N; nn++) {
				double sum1 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum1 += L.get(nn, k)*U.get(k, nn);
				}
				double ann = A.get(nn, n);
				if(Math.abs(ann-sum1) > 1e-12) {
					if(nn!=n) {
						int tmp = VP[nn];
						VP[nn] = VP[n];
						VP[n] = tmp;
						double vtmp;
						for(int c=1;c<=N;c++) {
							vtmp = A.get(n, c);
							A.set(n, c, A.get(nn, c));
							A.set(nn, c, vtmp);
						}
						for(int c=1;c<n;c++) {
							vtmp = L.get(n, c);
							L.set(n, c, L.get(nn, c));
							L.set(nn, c, vtmp);
						}
					}
					break;
				}
			}
			for(int j=n; j<=N; j++) {
				double sum1 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum1 += L.get(n, k)*U.get(k, j);
				}
				double anj = A.get(n, j);
				U.set(n, j, anj-sum1);
			}
			for(int i=n+1; i<=N; i++) {
				double sum2 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum2 += L.get(i, k)*U.get(k, n);
				}
				double ain = A.get(i, n);
				L.set(i, n, (ain-sum2)/U.get(n, n));
			}
			L.set(n, n, 1.0);
		}
		
		//P
		for(int i=1;i<=N;i++) {
			P.set(i, VP[i], 1.0);
		}
	}
	
	public static void LU(FullMatrix AA, FullMatrix L, FullMatrix U, SparseMatrix P) {
		FullMatrix A = AA.copy();
		int N = A.getRowDim();
		int[] VP = new int[N];
		for(int n=0; n<N; n++) {
			VP[n] = n;
		}
		double[][] dA = A.data;
		double[][] dL = L.data;
		double[][] dU = U.data;
		
		for(int n=0; n<N; n++) {
			//判断Unn的值是否为0，为0的话进行行对换，直到找到Unn不为0的行
			for(int nn=n; nn<N; nn++) {
				double sum1 = 0.0;
				for(int k=0; k<=n-1; k++) {
					sum1 += dL[nn][k]*dU[k][nn];
				}
				double ann = dA[nn][n];
				if(Math.abs(ann-sum1) > 1e-12) {
					if(nn!=n) {
						int tmp = VP[nn];
						VP[nn] = VP[n];
						VP[n] = tmp;
						double vtmp;
						for(int c=0;c<N;c++) {
							vtmp = dA[n][c];
							dA[n][c] = dA[nn][c];
							dA[nn][c] = vtmp;
						}
						for(int c=0;c<n;c++) {
							vtmp = dL[n][c];
							dL[n][c] = dL[nn][c];
							dL[nn][c] = vtmp;
						}
					}
					break;
				}
			}
			for(int j=n; j<N; j++) {
				double sum1 = 0.0;
				for(int k=0; k<=n-1; k++) {
					sum1 += dL[n][k]*dU[k][j];
				}
				double anj = dA[n][j];
				dU[n][j] = anj-sum1;
			}
			for(int i=n+1; i<N; i++) {
				double sum2 = 0.0;
				for(int k=0; k<=n-1; k++) {
					sum2 += dL[i][k]*dU[k][n];
				}
				double ain = dA[i][n];
				dL[i][n] = (ain-sum2)/dU[n][n];
			}
			dL[n][n] = 1.0;
		}
		
		//P
		for(int i=1;i<=N;i++) {
			P.set(i, VP[i-1]+1, 1.0);
		}
	}	
	/**
	 * U*x=f
	 * 
	 * @param U is an upper triangular matrix
	 * @param x
	 * @param f
	 */
	public static Vector solveUx(Matrix U, Vector x, Vector f) {
		int n = U.getRowDim();
		for(int i=n;i>0;i--) {
			double sum = 0.0;
			for(int j=n;j>i;j--) {
				sum += U.get(i, j)*x.get(j);
			}
			double xi = (f.get(i) - sum)/U.get(i, i);
			x.set(i, xi);
		}
		return x;
	}
	public static FullVector solveUx(FullMatrix U, FullVector x, FullVector f) {
		int N = U.getRowDim();
		double[][] dU = U.data;
		double[] dx = x.data;
		double[] df = f.data;
		for(int i=N-1;i>=0;i--) {
			double sum = 0.0;
			for(int j=N-1;j>i;j--) {
				sum += dU[i][j]*dx[j];
			}
			double xi = (df[i] - sum)/dU[i][i];
			dx[i] = xi;
		}
		return x;
	}	
	
	/**
	 * L*x = f
	 * 
	 * @param L is a lower triangular matrix with ones on its diagonal
	 * @param x
	 * @param f
	 * @return
	 */
	public static Vector solveLx(Matrix L, Vector x, Vector f) {
		int n = L.getRowDim();
		for(int i=1;i<=n;i++) {
			double sum = 0.0;
			for(int j=1;j<i;j++) {
				sum += L.get(i, j)*x.get(j);
			}
			double xi = f.get(i) - sum;
			x.set(i, xi);
		}
		return x;		
	}

	public static FullVector solveLx(FullMatrix L, FullVector x, FullVector f) {
		int N = L.getRowDim();
		double[][] dL = L.data;
		double[] dx = x.data;
		double[] df = f.data;
		for(int i=0;i<N;i++) {
			double sum = 0.0;
			for(int j=0;j<i;j++) {
				sum += dL[i][j]*dx[j];
			}
			double xi = df[i] - sum;
			dx[i] = xi;
		}
		return x;
	}
	
	/**
	 * P*x = f
	 * 
	 * @param P is a permutation matrix
	 * @param x
	 * @param f
	 * @return
	 */
	public static Vector solvePx(SparseMatrix P, Vector x, Vector f) {
		Map<Integer, Map<Integer, Double>> p = P.getAll();
		int n = P.getRowDim();
		for(int i=1;i<=n;i++) {
			Map<Integer, Double> row = p.get(i);
			for(Entry<Integer, Double> c : row.entrySet()) {
				x.set(c.getKey(),f.get(i));
			}
		}
		return x;
	}
	public static FullVector solvePx(SparseMatrix P, FullVector x, FullVector f) {
		Map<Integer, Map<Integer, Double>> p = P.getAll();
		int n = P.getRowDim();
		for(int i=1;i<=n;i++) {
			Map<Integer, Double> row = p.get(i);
			for(Entry<Integer, Double> c : row.entrySet()) {
				x.data[c.getKey()-1] = f.data[i-1];
			}
		}
		return x;
	}
	
	public static Vector solve(Matrix A, Matrix L, Matrix U, SparseMatrix P,
			Vector x, Vector f) {
		LU(A, L, U, P);
		solvePx(P,x,f);
		Vector x2 = x.copy();
		solveLx(L,x2,x);
		solveUx(U,x,x2);
		return x;
	}
	
	public static FullVector solve(FullMatrix A, FullMatrix L, FullMatrix U, SparseMatrix P,
			FullVector x, FullVector f) {
		LU(A, L, U, P);
		solvePx(P,x,f);
		FullVector x2 = x.copy();
		solveLx(L,x2,x);
		solveUx(U,x,x2);
		return x;
	}
	
	
	public static void test1() {
		SparseMatrix A = new SparseMatrix(3,3);
		SparseMatrix L = new SparseMatrix(3,3);
		SparseMatrix U = new SparseMatrix(3,3);
		SparseMatrix P = new SparseMatrix(3,3);
		
		double[][] data = {{8,2,9},{4,9,4},{6,7,9}};
		//double[][] data = {{1,0,0},{0,0,2},{0,1,-1}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		LU(A,L,U,P);
		A.print();
		L.print();
		U.print();
		P.print();

		SparseVector x = new SparseVector(3);
		SparseVector f = new SparseVector(3,1.0);
		solve(A,L,U,P,x,f);
		x.print();
		
		SparseVector y = new SparseVector(3);
		A.mult(x, y);
		y.print();		
	}
	
	public static void test2() {
		SparseMatrix A = new SparseMatrix(3,3);
		SparseMatrix L = new SparseMatrix(3,3);
		SparseMatrix U = new SparseMatrix(3,3);
		SparseMatrix P = new SparseMatrix(3,3);
		
		double[][] data = {{8,2,9},{4,9,4},{6,7,9}};
		//double[][] data = {{1,0,0},{0,0,2},{0,1,-1}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		
		FullMatrix fA = new FullMatrix(A);
		FullMatrix fL = new FullMatrix(L);
		FullMatrix fU = new FullMatrix(U);
		LU(fA,fL,fU,P);
		fA.print();
		fL.print();
		fU.print();
		P.print();

		FullVector x = new FullVector(3);
		FullVector f = new FullVector(3,1.0);
		solve(fA,fL,fU,P,x,f);
		x.print();
		
		SparseVector y = new SparseVector(3);
		A.mult(x.getSparseVector(), y);
		y.print();
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		test1();
		System.out.println("-----------------------------");
		test2();
	}

}
