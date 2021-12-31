package org.cnrs.crbm.lib.math;
import org.apache.commons.lang3.ArrayUtils;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 11/19/2014.
 */
public class PeriodFinder {


    public double findPeriodicity(double[] input) {
        double period = 0.0;

        //for sort
        int[] peaks = this.getBottomsWithSize(input, 5);


//        double deff = 0.0;
//        for (int i = 0; i < peakLocations.length; i++) {
//            //if (input[peakLocations[i]] < 6.0)
//            //storePeaks.add(peakLocations[i]);
//            // System.out.println(peakLocations[i]);
//        }
//        //
//        double[] y = new double[peakLocations.length];
//        for (int i = 0; i < peakLocations.length; i++) {
//            y[i] = input[peakLocations[i]];
//
//            // System.out.print(NumberFormatUtils.format(y[i]) + " ");
//
//        }
//
//
//        // find long period
//        int[] peaks = this.getBottomsWithSize(y, 2);
        List<Integer> storePeaks = new ArrayList<Integer>();
        for (int i = 0; i < peaks.length; i++) {
            // System.out.print(peaks[i] + " ");

            if (input[peaks[i]] < 5.0)
                storePeaks.add(peaks[i]);


            // System.out.print(peakLocations[peaks[i]] + " ");
        }
        // System.out.println();



        peaks = ArrayUtils.toPrimitive((Integer[]) storePeaks.toArray(new Integer[storePeaks.size()]));
        double deff = 0.0;
        for (int i = 0; i < peaks.length - 1; i++) {

            //System.out.println(peakLocations[peaks[i + 1]] - peakLocations[peaks[i]]);
            deff += peaks[i + 1] - peaks[i];

        }

        double longLen = 0.0;
        if (peaks.length >= 2)
            longLen = deff / (peaks.length - 1);


        period = longLen;

        return period;
    }

    public int[] getBottomsWithSize(double[] signals, int windowsize) {
        // convert signals to peaks signals
        double[] signals_down = new double[signals.length];
        for (int i = 0; i < signals.length; i++) {
            signals_down[i] = 10 - signals[i];
            //System.out.print(NumberFormatUtils.format(signals_down[i]) + ",");
        }
        PeakDetector detector = new PeakDetector(signals_down);
        //windowsize = 3;
        int stringency = 1;

        int[] peakLocations = detector.process(windowsize, stringency);
        return peakLocations;
    }
}
