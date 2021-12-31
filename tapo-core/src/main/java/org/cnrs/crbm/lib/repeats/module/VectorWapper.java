package org.cnrs.crbm.lib.repeats.module;

import org.apache.commons.math3.ml.clustering.Clusterable;

public class VectorWapper implements Clusterable {
	private double[] points;
	private ProVector vector;
	private int index;

	public VectorWapper(ProVector vector, double angle, int index) {
		this.vector = vector;
		this.index = index;
		this.points = new double[] { angle, 0 };
	}

	public ProVector getVector() {
		return vector;
	}

	public double[] getPoint() {
		return points;
	}

	public int getIndex() {
		return index;
	}

}