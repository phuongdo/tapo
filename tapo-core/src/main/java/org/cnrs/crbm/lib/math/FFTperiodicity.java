package org.cnrs.crbm.lib.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.jtransforms.fft.DoubleFFT_1D;

public class FFTperiodicity {

    public static void main(String[] args) throws Exception {
        FFTperiodicity fft = new FFTperiodicity();

        double[] inputStart = new double[300];
        for (int x = 0; x < inputStart.length; x++) {
            inputStart[x] = Math.sin(x);
        }

        inputStart = Filter.filter(inputStart);

        //System.out.println(inputStart.length);
        //fft.show(inputStart);


        int winsize = 40;
        //fft.findPeriodicityWithWindow(inputStart, winsize);

        System.out.println(fft.findPeriodicity(inputStart));

    }

    public double findPeak(double[] inputStart) {

        PeakDetector detector = new PeakDetector(inputStart);
        DescriptiveStatistics stats = new DescriptiveStatistics();
        // int repeatLeng = 0;
        int[] peakLocations = detector.process(5, 1);
        double[] peakTime = new double[peakLocations.length];
        double[] peakVal = new double[peakLocations.length];
        if (peakLocations.length >= 3) {
            for (int i = 0; i < peakLocations.length; i++) {
                peakTime[i] = peakLocations[i];
                peakVal[i] = inputStart[peakLocations[i]];
                // System.out
                // .println("At t=" + peakTime[i] + " - Val = " + peakVal[i]);

                stats.addValue(peakVal[i]);
            }

            double P_avg = stats.getMean();
            double P_sd = stats.getStandardDeviation();
            stats = new DescriptiveStatistics();
            for (int i = 0; i < peakLocations.length; i++) {
                peakTime[i] = peakLocations[i];
                peakVal[i] = inputStart[peakLocations[i]];

                if (P_avg - P_sd / 2 <= peakVal[i]
                        && peakVal[i] <= P_avg + P_sd / 2)
                    stats.addValue(peakVal[i]);
                // System.out
                // .println("At t=" + peakTime[i] + " - Val = " + peakVal[i]);

//                stats.addValue(peakVal[i]);
            }
            return stats.getMean();
        }
        return 0;

    }

    public List<Double> findPeriodicityWithConstantSignal(double[] inputStart,
                                                          int winsize) {

        // Get a DescriptiveStatistics instance
        DescriptiveStatistics stats = new DescriptiveStatistics();
        List<Double> posibleLengs = new ArrayList<Double>();

        // List<Double> listPeriods = new ArrayList<Double>();
        for (int i = 0; i < inputStart.length - winsize; i = i + winsize) {
            double[] v = Arrays.copyOfRange(inputStart, i, i + winsize);
            double period = this.findPeak(v);
            // listPeriods.add(period);
            if (period > 10)
                posibleLengs.add(period);

        }

        // double P_avg = stats.getMean();
        // double P_sd = stats.getStandardDeviation();
        // stats = new DescriptiveStatistics();
        //
        // for (Double period : listPeriods) {
        // if (P_avg - P_sd / 2 <= period && period <= P_avg + P_sd / 2)
        // posibleLengs.add(period);
        // }

        return posibleLengs;

    }

    public List<Double> findPeriodicityWithWindow(double[] inputStart,
                                                  int winsize) {

        // Get a DescriptiveStatistics instance
        DescriptiveStatistics stats = new DescriptiveStatistics();
        List<Double> posibleLengs = new ArrayList<Double>();

        List<Double> listPeriods = new ArrayList<Double>();
        for (int i = 0; i < inputStart.length - 5 && i + winsize <= inputStart.length; i = i + 5) {
            double[] v = Arrays.copyOfRange(inputStart, i, i + winsize);
            double period = this.findPeriodicity(v);
            // System.out.println(this.findPeak(v));
            // stats.addValue(period);
            listPeriods.add(period);
            // posibleLengs.add(this.findPeak(v));


        }

        Double[] periods = listPeriods.toArray(new Double[listPeriods.size()]);

        int[] labels = new int[periods.length];
        for (int i = 0; i < labels.length; i++)
            labels[i] = -1;

        int label = 0;
        for (int i = 0; i < periods.length; i++) {

            if (labels[i] == -1) {
                Double periodRf = periods[i];
                labels[i] = label;

                for (int j = i + 1; j < periods.length; j++) {
                    if (labels[j] == -1)
                        if (periodRf - 5 <= periods[j]
                                && periods[j] <= periodRf + 5) {
                            labels[j] = label;
                        }
                }

                label++;
            }

        }
        Map<Integer, Double> C = this.getAvgLeng(labels, periods);

        for (Entry<Integer, Double> entry : C.entrySet()) {
            posibleLengs.add(entry.getValue());

        }
        // double P_avg = stats.getMean();
        // double P_sd = stats.getStandardDeviation();
        // stats = new DescriptiveStatistics();
        //
        // for (Double period : periods) {
        // if (P_avg - P_sd / 2 <= period && period <= P_avg + P_sd / 2)
        // posibleLengs.add(period);
        // }
        return posibleLengs;

    }

    public Map<Integer, Double> getAvgLeng(int[] labels, Double[] periods) {
        Map<Integer, Double> avg = new HashMap<Integer, Double>();
        Map<Integer, Integer> count = new HashMap<Integer, Integer>();
        for (int i = 0; i < labels.length; i++) {

            int label = labels[i];
            if (count.containsKey(label)) {

                count.put(label, count.get(label) + 1);
            } else {
                count.put(label, 1);

            }

        }
        for (int i = 0; i < labels.length; i++) {

            int label = labels[i];
            if (avg.containsKey(label)) {

                avg.put(label, avg.get(label) + periods[i]);
            } else {
                avg.put(label, periods[i]);

            }

        }

        for (Entry<Integer, Double> entry : avg.entrySet()) {
            int label = entry.getKey();

            if (count.containsKey(label)) {

                avg.put(label, (double) avg.get(label) / count.get(label));
            }

        }

        return avg;
    }


    public double findPeriodicity(double[] inputStart) {

        double[] complex = this.transformNew(inputStart);
        int n = inputStart.length;
        double[] power = new double[n];


        for (int i = 0; i < n; i++) {

            double re = complex[2 * i];
            double im = complex[2 * i + 1];
            power[i] = (re * re + im * im) / Math.sqrt(n);

        }

        double nyquist = 0.5;
        int mid_n = n / 2;

        // (1:n/2)/(n/2)*nyquist;
        double[] freq = new double[mid_n];
        for (int i = 0; i < mid_n; i++) {
            freq[i] = ((double) i / mid_n) * nyquist;
        }
        double[] period = new double[mid_n];
        for (int i = 0; i < mid_n; i++) {
            period[i] = (double) 1 / freq[i];
        }
        power = Arrays.copyOfRange(power, 1, mid_n + 1);
        period = Arrays.copyOfRange(period, 1, mid_n + 1);
//        show(ts);
//        System.out.println();
//        show(period);

        int index = this.findMaxIndexArg(power, period);
        return period[index];

    }

    public int findMaxIndexArg(double[] v, double[] p) {
        double max = Double.MIN_VALUE;

        int index = 0;
        for (int i = v.length - 1; i > 0 && p[i] < 100; i--) {

            if (max < v[i]) {
                max = v[i];
                index = i;
            }
        }
        //System.err.println(max);
        return index;

    }

    public void show(double[] v) {

        for (int i = 0; i < v.length; i++)
            System.out.print(NumberFormatUtils.format(v[i]) + ", ");
        //System.out.print(v[i]) + " , ");
    }


    public double[] transformNew(double[] input) {


        DoubleFFT_1D fftDo = new DoubleFFT_1D(input.length);
        double[] fft = new double[input.length * 2];
        System.arraycopy(input, 0, fft, 0, input.length);
        fftDo.realForwardFull(fft);
        return fft;
    }

    public double[] transform(double[] inputStart) {

        // copy array
        int length = (int) Math.pow(2,
                Math.ceil(Math.log(inputStart.length) / Math.log(2)));
        double[] input = new double[length];
        for (int i = 0; i < inputStart.length; i++) {
            input[i] = inputStart[i];
        }
        double[] tempConversion = new double[input.length];

        FastFourierTransformer transformer = new FastFourierTransformer(
                DftNormalization.STANDARD);
        try {
            Complex[] complx = transformer.transform(input,
                    TransformType.FORWARD);

            for (int i = 0; i < complx.length; i++) {
                double rr = (complx[i].getReal());
                double ri = (complx[i].getImaginary());
                tempConversion[i] = Math.sqrt((rr * rr) + (ri * ri));
            }
        } catch (IllegalArgumentException e) {
            System.out.println(e);
        }

        return tempConversion;
    }
}
