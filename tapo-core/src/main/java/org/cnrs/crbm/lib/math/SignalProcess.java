package org.cnrs.crbm.lib.math;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.cnrs.crbm.lib.analysis.PowerRepeatVsNoRepeat;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by pdoviet on 11/3/2014.
 */
public class SignalProcess {


    public SignalPeriod costanceDetect(SignalPeriod signal) {

        /* periodic = -1 means that this signals is not constance signal */
        double periodic = -1;


        double[] s = signal.getS();
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = 0; i < s.length; i++) {

            stats.addValue(s[i]);
        }


        double P_avg = stats.getMean();
        double P_sd = stats.getStandardDeviation();
        stats = new DescriptiveStatistics();
        int match = 0;
        for (int i = 0; i < s.length; i++) {
            if (P_avg - P_sd / 2 <= s[i]
                    && s[i] <= P_avg + P_sd / 2) {
                stats.addValue(s[i]);
                match++;
            }


        }

        if ((double) match / s.length > 0.8)
            periodic = stats.getMean();
        else
            periodic = -1;
        //return stats.getMean();

        signal.setPeriod(periodic);
        return signal;
    }

    public SignalPeriod periodDetect(SignalPeriod signal) throws Exception {

        /* periodic = -1 means that this signals is not constance signal */
        double periodic = -1;
        FFTperiodicity fft = new FFTperiodicity();
        double[] smootSignal = Filter.filter(signal.getS());
        periodic = fft.findPeriodicity(smootSignal);
        signal.setPeriod(periodic);
        return signal;
    }


    public double predictPeriod(double[] his) throws Exception {
        FFTperiodicity fft = new FFTperiodicity();
        double periodic = fft.findPeriodicity(his);
        return periodic;
    }

    PowerRepeatVsNoRepeat powerRepeatVsNoRepeat = new PowerRepeatVsNoRepeat();

    public List<SignalPeriod> predictPeriodicties(double[] his) throws Exception {
        SignalProcess signalProcess = new SignalProcess();
        List<SignalPeriod> posibleLengs = new ArrayList<SignalPeriod>();
        int[] wlist = new int[]{200};
//        if (his.length > 500)
//            wlist = new int[]{his.length / 3};
        for (Integer winsize : wlist) {
            for (int i = 0; i < his.length; ) {
                try {
                    int start = i;
                    int end = i + winsize;

                    if (end > his.length)
                        end = his.length;
                    double[] v = Arrays.copyOfRange(his, start, end);

                    if ((double) (end - start + 1) / winsize < 0.2)
                        break;

                    //System.out.println("XXX");
                    SignalPeriod signalPeriod = new SignalPeriod(start, end, v);
                    // process v here
                    signalPeriod = signalProcess.costanceDetect(signalPeriod);
                    double period = -1;
                    if (signalPeriod.getPeriod() == -1) {
                        signalPeriod = signalProcess.periodDetect(signalPeriod);
                    }

                    double sgScore = this.correlationWithLag(v, (int) signalPeriod.getPeriod());
                    signalPeriod.setSgScore(sgScore);
                    posibleLengs.add(signalPeriod);

                    i = i + winsize / 2;
                } catch (Exception ex) {
                    ///
//                    ex.printStackTrace();
                    //System.out.println(ex);
                    i = i + winsize / 2;
                }

            }
        }

        return posibleLengs;
    }

    public double correlationWithLag(double[] signal, int lag) {

        double corr = 0.0;
        int signallenth = signal.length;
        SpearmansCorrelation correlation = new SpearmansCorrelation();
        //PearsonsCorrelation correlation = new PearsonsCorrelation();
        double[] xArray = new double[signallenth - lag];
        double[] yArray = new double[signallenth - lag];

        for (int i = 0; i < signallenth - lag; i++) {
            xArray[i] = signal[i];
            yArray[i] = signal[i + lag];
        }
        corr = correlation.correlation(xArray, yArray);
        return corr;
    }



//    public static void main(String[] args) {
//
//
//        double[] his = new double[]{0, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 100, 23, 32, 23};
//        SignalProcess signalProcess = new SignalProcess();
//        System.out.println(signalProcess.costanceDetect(his));
//
//
//
//    }


}
