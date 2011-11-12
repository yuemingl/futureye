package edu.uta.futureye.algebra;

import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;

public class FullMatrix implements AlgebraMatrix {
	protected double[][] data = null;
	protected int rowDim;
	protected int colDim;

	public FullMatrix() {
		
	}
	
	public FullMatrix(int nRow, int nCol) {
		this.rowDim = nRow;
		this.colDim = nCol;
		
		this.data = new double[nRow][];
		for(int r=0; r<nRow; r++) {
			this.data[r] = new double[nCol];
			for(int i=0;i<nCol;i++) {
				data[r][i] = 0.0;
			}
		}
	}
	
	public FullMatrix(SparseMatrix sMat) {
		this.rowDim = sMat.getRowDim();
		this.colDim = sMat.getColDim();
		
		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
		this.data = new double[this.rowDim][];
		for(int r=0; r<this.rowDim; r++) {
			this.data[r] = new double[this.colDim];
			for(int i=0;i<this.colDim;i++) {
				data[r][i] = 0.0;
			}
			Map<Integer, Double> row = m.get(r+1);
			if(row != null) {
				for(Entry<Integer,Double> e : row.entrySet()) {
					this.data[r][e.getKey()-1] = e.getValue();
				}
			}
		}
	}
	
	@Override
	public int getRowDim() {
		return rowDim;
	}

	@Override
	public int getColDim() {
		return colDim;
	}

	@Override
	public void mult(AlgebraVector x, AlgebraVector y) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mult(AlgebraMatrix B, AlgebraMatrix C) {
		// TODO Auto-generated method stub

	}

	@Override
	public AlgebraMatrix getTrans() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void print() {
		int nRow = this.rowDim;
		int nCol = this.colDim;
		for(int i=0;i<nRow;i++) {
			for(int j=0;j<nCol;j++) {
				System.out.print(String.format("%8.6f   ", this.data[i][j]));
			}
			System.out.println();
		}
		System.out.println();

	}
	
	FullMatrix copy() {
		FullMatrix rlt = new FullMatrix();
		int nRow = this.rowDim;
		int nCol = this.colDim;
		rlt.rowDim = nRow;
		rlt.colDim = nCol;
		rlt.data = new double[nRow][];
		for(int i=0;i<nRow;i++) {
			rlt.data[i] = new double[nCol];
			for(int j=0;j<nCol;j++) {
				rlt.data[i][j] = this.data[i][j];
			}
		}
		return rlt;
	}

}
