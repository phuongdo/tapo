package org.cnrs.crbm.lib.repeats.module;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Created by pdoviet on 4/24/2015.
 */
public class Dihdral {
    public static double calcDihedral(Vector3D p1, Vector3D p2, Vector3D p3, Vector3D p4) {

        // p1-p2-p3-p4 and we need dihedral around the p2-p3 bond
        // this is the angle between two planes: (p1 p2 p3) and (p2 p3 p4)
        Vector3D v1 = p1.subtract(p2); // points from p2 to p1
        Vector3D v2 = p3.subtract(p2); // points from p2 to p3
        Vector3D n1 = v1.crossProduct(v2); // perpendicular to (p1 p2 p3) plane
        Vector3D v3 = p2.subtract(p3); // points from p3 to p2
        Vector3D v4 = p4.subtract(p3); // points from p3 to p4
        Vector3D n2 = v3.crossProduct(v4); // perpendicular to (p2 p3 p4) plane

        //double da = n1.angle(n2); // angle between two vectors perpendicular

        double da = Math.toDegrees(Vector3D.angle(n1, n2));

        // to (p1 p2 p3) and (p2 p3 p4) planes,
        // respectively.
        // angle between the two vectors is equal
        // to the angle between the two planes

        return Vector3D.dotProduct(n1, v4) < 0 ? -da : da;

    }


}
