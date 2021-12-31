package org.cnrs.crbm.lib.math;

import org.cnrs.crbm.lib.repeats.Point;

/**
 * Created by pdoviet on 2/5/2015.
 */
public class Function {

    // y = ax+b;
    double a;
    double b;

    public double getA() {
        return a;
    }

    public void setA(double a) {
        this.a = a;
    }

    public double getB() {
        return b;
    }

    public void setB(double b) {
        this.b = b;
    }

    @Override
    public String toString() {
        return "y = " + a + "x+" + b;
    }


    public double getValue(Point x) {

        return a * x.getX() + b;

    }
}
