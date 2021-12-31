package org.cnrs.crbm.lib.math;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

public class VectorCals {
	/**
	 * Projects this original vector onto another vector
	 * 
	 * @param other
	 *            The vector to project this vector onto
	 * @return A new vector with the projection result
	 */
	public static Vector3D project(Vector3D orignal, Vector3D other) {
		double n = orignal.dotProduct(other); // A . B
		double d = other.getNormSq(); // |B|^2
		return new Vector3D(other.getX(), other.getY(), other.getZ())
				.scalarMultiply(n / d);
	}
}
