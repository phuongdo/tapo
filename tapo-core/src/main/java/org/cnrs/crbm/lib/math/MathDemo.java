package org.cnrs.crbm.lib.math;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.geometry.euclidean.twod.Line;

public class MathDemo {

	public static void main(String[] args) {

		Vector3D vector3d = new Vector3D(1, 2, 3);
		Vector3D vectorSu = vector3d.add(new Vector3D(2, 1, 1));
		System.out.println(vectorSu);
		

	}

}
