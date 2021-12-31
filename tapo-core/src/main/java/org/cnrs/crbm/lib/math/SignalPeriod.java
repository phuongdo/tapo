package org.cnrs.crbm.lib.math;

import org.cnrs.crbm.lib.utils.NumberFormatUtils;

/**
 * Created by pdoviet on 11/4/2014.
 */
public class SignalPeriod {


    private int start;
    private int end;
    private double period;
    private  double sgScore;

    public double getSgScore() {
        return sgScore;
    }

    public void setSgScore(double sgScore) {
        this.sgScore = sgScore;
    }
    public double[] getS() {
        return s;
    }

    public double getPeriod() {
        return period;
    }

    public void setPeriod(double period) {
        this.period = period;
    }

    public void setS(double[] s) {
        this.s = s;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    private double[] s;

    public SignalPeriod(int start, int end, double[] s) {
        this.start = start;
        this.end = end;
        this.s = s;
    }
    public SignalPeriod(){}

    @Override
    public String toString() {
        return start + "-" + end + "\t" + NumberFormatUtils.format(period);
    }
}
