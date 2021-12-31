package org.cnrs.crbm.maths;

import org.apache.commons.math3.ml.clustering.Clusterable;

/**
 * Created by pdoviet on 8/17/2015.
 */
public class LocationWrapper implements Clusterable {
    private double[] points;


    public LocationWrapper(double avgLen) {

        this.points = new double[]{avgLen, 0};
    }


    public double[] getPoint() {
        return points;
    }
}