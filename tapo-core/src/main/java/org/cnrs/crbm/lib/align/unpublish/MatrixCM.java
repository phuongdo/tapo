package org.cnrs.crbm.lib.align.unpublish;

public class MatrixCM {

	private int[][] matrix;
	private int size;

	public MatrixCM(int size) {
		this.setSize(size);
		this.matrix = new int[size][size];
	}

	public int getCell(int r, int c) {
		return matrix[r][c];
	}

	public void setCell(int r, int c, int val) {
		matrix[r][c] = val;
	}

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public MatrixCM getSubMatrixCM(int start, int end) {

		MatrixCM sub = new MatrixCM(end - start + 1);
		for (int i = start; i <= end; i++) {

			for (int j = i + 1; j <= end; j++) {
				sub.setCell(i - start, j - start, this.getCell(i, j));
			}
		}

		return sub;
	}

	public boolean isOverllap(MatrixCM matrixCM) {
		if (matrixCM.getSize() != this.size)
			return false;
		// check overllap
		for (int i = 0; i < this.size; i++) {
			for (int j = i + 1; j < this.size; j++) {

				if ((this.getCell(i, j) != matrixCM.getCell(i, j)))
					return false;
			}
		}
		return true;

	}
}
