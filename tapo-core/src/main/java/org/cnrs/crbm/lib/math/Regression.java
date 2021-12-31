package org.cnrs.crbm.lib.math;

import org.cnrs.crbm.lib.repeats.Point;

/**
 * Created by pdoviet on 2/5/2015.
 */
public class Regression {

    public static void main(String[] args) {


        Point a = new Point();
        a.setX(3);
        a.setY(7);

        Point b = new Point();
        b.setX(5);
        b.setY(11);

        Function f = new Regression().getFunction(a,b);
        System.out.println(f);



    }


    public static Function getFunction(Point x1, Point x2) {

        //y = mx+b;
        double a = (x2.getY() - x1.getY()) / (x2.getX() - x1.getX());
        double b = x1.getY() - a * x1.getX();
        Function function = new Function();
        function.setA(a);
        function.setB(b);
        return function;

    }


}
