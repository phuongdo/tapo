package org.cnrs.crbm.lib.analysis;

import com.google.common.math.IntMath;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.math.FFTperiodicity;
import org.cnrs.crbm.lib.math.Filter;
import org.cnrs.crbm.lib.math.PeakDetector;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by pdoviet on 3/9/2015.
 */
public class PowerRepeatVsNoRepeat {

    public static void main(String[] args) throws Exception {

        new PowerRepeatVsNoRepeat().run();

//        new PowerRepeatVsNoRepeat().analyse("1ap7", "A", "0");
    }

    Superimposer superimposer = new Superimposer();

    public void analyse(String code, String chain, String type) {

        try {
            RepeatFinder finder = new RepeatFinder(code, chain);
            Atom[] atoms = finder.getFeatures().getAtoms();
            //RMSD signals
            String strSS = finder.getFeatures().getStrSS();

            //MutilAlign pairAlign = new MutilAlign();


//            this.show(his);
//            System.out.println();
//            System.exit(1);


            //System.out.println(predictPeriod);

            //fft.show(his);
            //fft.show(period);
            //System.out.println();
            //fft.show(power);
            //System.out.println();
            //System.out.println(period[index]);


            SignalData signalData = this.processStructure(atoms);
            System.out.println(code + "_" + chain + "\t" + type + "\t" + atoms.length + "\t" + NumberFormatUtils.format(signalData.zScore) + "\t" + NumberFormatUtils.format(signalData.phi) + "\t" + NumberFormatUtils.format(signalData.variation) + "\t" + NumberFormatUtils.format(signalData.corr));
        } catch (Exception ex) {
            //ex.printStackTrace();
        }

    }


    public double getSignalScore(Features features) throws Exception {
        double signalScore = 0.0;
        Atom[] atoms = features.getAtoms();
        SignalData signalDataRMSD = this.processStructure(atoms);
        SignalData signalDataCM = this.proCessCMSignal(features);
        signalScore = Math.max(signalDataRMSD.corr,signalDataCM.corr);
        return signalScore;
    }


    public  SignalData proCessCMSignal(Features features) throws Exception {
        double[] his = features.getCmHistos();
        double[] smootSignal = Filter.filter(his);
        double maxZ = Double.MIN_VALUE;
        double maxV = 0.0;
        double var = 0.0;
        double phi = 0.0;
        double maxCorr = 0.0;
        int winsize = 100;
        for (int i = 0; i < smootSignal.length - winsize; i = i + winsize / 4) {
            try {
                int start = i;
                int end = i + winsize;

                if (end > smootSignal.length)
                    end = smootSignal.length;

                double[] v = Arrays.copyOfRange(smootSignal, start, end);

                if (v.length < winsize)
                    continue;
                SignalData signalData = this.analysisSignal(v);
                //System.out.println(signalData.zScore);

                if (maxCorr < signalData.corr) {
                    maxCorr = signalData.corr;
                    maxV = signalData.variation;
                    maxZ = signalData.zScore;
                    var = signalData.variation;
                    phi = signalData.phi;
                }

            } catch (Exception ex) {
                //
            }

        }
        SignalData signalData = new SignalData();
        signalData.corr = maxCorr;
        signalData.variation = var;
        signalData.zScore = maxZ;
        signalData.phi = phi;
        return signalData;

    }

    public SignalData analysisSegmentAtoms(Atom[] atoms) {
        List signals = new ArrayList<Double>();
        int winsize = 20;
        int t = 10;
        Atom[] seedAtoms = Fragement.getFragementsofAtoms(atoms, t, t
                + winsize);
        for (int j = 0; j < atoms.length - winsize; j++) {
            Atom[] compareAtoms = Fragement.getFragementsofAtoms(atoms, j,
                    j + winsize);

            double rmsd = 10;
            try {
                rmsd = superimposer.superimposeSimple(seedAtoms,
                        compareAtoms);
                signals.add(rmsd);
            } catch (StructureException e) {
                e.printStackTrace();
            }

        }

        double[] his = ArrayUtils.toPrimitive((Double[]) signals.toArray(new Double[signals.size()]));
//            double[] his = finder.getCmHistos();
        try {
            his = Filter.filter(his);
        } catch (Exception e) {
            e.printStackTrace();
        }
        //his = VectorSmoother.smoothWithDegreeFive(his);
        //System.out.println(this.correlationWithLag(his, 21));
        //this.show(his);
        SignalData signalData = this.analysisSignal(his);
        return signalData;
    }

    public SignalData processStructure(Atom[] atoms) {

        double maxZ = Double.MIN_VALUE;
        double maxV = 0.0;
        double var = 0.0;
        double phi = 0.0;
        double maxCorr = 0.0;
        int winsize = 100;
        for (int i = 0; i < atoms.length - winsize; i = i + winsize / 4) {
            try {
                int start = i;
                int end = i + winsize;

                if (end > atoms.length)
                    end = atoms.length;

                Atom[] v = Arrays.copyOfRange(atoms, start, end);

                if (v.length < winsize)
                    continue;
                SignalData signalData = this.analysisSegmentAtoms(v);
                //System.out.println(signalData.zScore);

                if (maxCorr < signalData.corr) {
                    maxCorr = signalData.corr;
                    maxV = signalData.variation;
                    maxZ = signalData.zScore;
                    var = signalData.variation;
                    phi = signalData.phi;

                }

            } catch (Exception ex) {
                //
            }

        }

        SignalData signalData = new SignalData();
        signalData.corr = maxCorr;
        signalData.variation = var;
        signalData.zScore = maxZ;
        signalData.phi = phi;
        return signalData;


    }


    private class SignalData {

        public double zScore = 0.0;
        public double phi = 0.0;
        public double variation = 0.0;
        public double corr = 0.0;


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


    public SignalData analysisSignal(double[] his) {
        SignalData signalData = new SignalData();
        FFTperiodicity fft = new FFTperiodicity();
        int n = his.length;
        int mid_n = IntMath.divide(n, 2, RoundingMode.UP);

        // mid_n = n;
        double[] complex = fft.transformNew(his);
        double[] power = new double[n];


        for (int i = 0; i < n; i++) {

            double re = complex[2 * i];
            double im = complex[2 * i + 1];
            power[i] = (re * re + im * im) / (n);

        }
        double nyquist = 0.5;


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

//            period = Arrays.copyOfRange(period, 1, mid_n);


        /**
         * PHI = [1,...,5]
         */
        double PHI = 3;

        DescriptiveStatistics stats = new DescriptiveStatistics();

        for (int i = 0; i < power.length; i++) {
            stats.addValue(power[i]);
        }

        double P_avg = stats.getMean();
        double P_sd = stats.getStandardDeviation();
        double[] z_score = new double[mid_n];

        for (int i = 0; i < power.length; i++) {
            stats.addValue(power[i]);
        }

        double N_PHI = 0;
        for (int i = 0; i < z_score.length; i++) {
            z_score[i] = (power[i] - P_avg) / P_sd;
            if (z_score[i] > PHI) {
                N_PHI = N_PHI + 1.0;
            }
        }

        double max = Double.MIN_VALUE;
        int index = 0;
        for (int i = z_score.length - 1; i > 0; i--) {

            if (max < z_score[i]) {
                max = z_score[i];
                index = i;
            }
        }

        double P_PHI = 100 * N_PHI / (mid_n);

        int predictPeriod = (int) period[index];

        // verify period


        PeakDetector detector = new PeakDetector(his);
        // get peaks

        int[] peaks = detector.process(3, 0.5);
        // get periods


        int N = 0;
        double D0 = 1;
        double delta = 0;
//        for (int i = 0; i < peaks.length - 1; i++) {
//
//            double d = his[peaks[i + 1]] - his[peaks[i]];
//            delta += 1 / (1 + Math.pow(predictPeriod - d / D0, 2));
//            N++;
//
//
//        }
        //System.out.println(predictPeriod);
        for (int i = 0; i < his.length; i++) {
            if (i < his.length && i + predictPeriod < his.length) {
                //System.out.println(his[i] - his[i + predictPeriod - 1]);
                delta += 1 / (1 + Math.pow((his[i] - his[i + predictPeriod - 1]) / D0, 2));
                N++;
            }
        }

        double variation = delta / N;
        signalData.corr = Math.abs(this.correlationWithLag(his, predictPeriod));
        signalData.variation = variation;
        signalData.zScore = max;
        signalData.phi = P_PHI;
        return signalData;

    }

    public SignalData processSignal(double[] signal) {

//        PeakDetector detector = new PeakDetector(signal);
//        // get peaks
        double maxZ = Double.MIN_VALUE;
        double maxV = 0.0;
        double var = 0.0;
        double phi = 0.0;
        double maxCorr = 0.0;
//        int[] peaks = detector.process(3, 0.5);
//        // get periods
//
//        int[] profiles = new int[signal.length];
//        if (peaks.length > 1) {
//            int[] periods = new int[peaks.length - 1];
//            for (int i = 1; i < peaks.length - 1; i++) {
//
//                periods[i] = peaks[i + 1] - peaks[i];
//                // profiles[peaks[i]] = periods[i];
//                for (int j = peaks[i]; j < peaks[i + 1]; j++) {
//                    profiles[j] = periods[i];
//                }
//            }
//
//
//            int[] labels = new int[periods.length];
//            for (int i = 0; i < labels.length; i++)
//                labels[i] = -1;
//
//            int label = 0;
//            for (int i = 0; i < periods.length; i++) {
//
//                if (labels[i] == -1) {
//                    int periodRf = periods[i];
//                    labels[i] = label;
//
//                    for (int j = i + 1; j < periods.length; j++) {
//                        if (labels[j] == -1)
//                            if (periodRf - 5 <= periods[j]
//                                    && periods[j] <= periodRf + 5) {
//                                labels[j] = label;
//                            }
//                    }
//
//                    label++;
//                }
//
//            }
//
//
////            double minVar = Double.MIN_VALUE;
////            double ma = Double.MIN_VALUE;
//
//            // extract the signals
//            for (int i = 0; i < label; i++) {
//
//                //System.out.println(i);
//                int start = 0;
//                int end = 0;
//                List<List<Integer>> regions = new ArrayList<List<Integer>>();
//                // extract all regions
//                for (int j = 0; j < labels.length; j++) {
//                    if (labels[j] == i) {
//                        start = j;
//                        while (j < labels.length
//                                && labels[j] == i) {
//                            j++;
//                        }
//                        j = end = j - 1;
//
//                        List<Integer> region = new ArrayList<Integer>();
//                        region.add(peaks[start]);
//                        region.add(peaks[end]);
//                        regions.add(region);
//                    }
//
//                }
//
//
//                for (List<Integer> region : regions) {
//                    int startR = region.get(0);
//                    int endR = region.get(region.size() - 1);
//                    if (endR - startR + 1 > 20) {
//                        //System.out.println(startR + ":" + endR);
//
//                        // this.show(Arrays.copyOfRange(signal, startR, endR + 1));
//                        SignalData signalData = this.analysisSignal(Arrays.copyOfRange(signal, startR, endR + 1));
//                        //System.out.println(signalData.zScore);
//
//                        if (maxZ < signalData.zScore) {
//                            maxZ = signalData.zScore;
//                            var = signalData.variation;
//                            phi = signalData.phi;
//
//                        }
//
//
//                    }
//                }
//
//            }
//
////            this.show(labels);
////            this.show(peaks);
//
//
//        }


        int winsize = 80;
        for (int i = 0; i < signal.length; i = i + winsize / 4) {

            try {
                int start = i;
                int end = i + winsize;

                if (end > signal.length)
                    end = signal.length;
                double[] v = Arrays.copyOfRange(signal, start, end);

                SignalData signalData = this.analysisSignal(v);
                //System.out.println(signalData.zScore);

                if (maxCorr < signalData.corr) {
                    maxCorr = signalData.corr;
                    maxV = signalData.variation;
                    maxZ = signalData.zScore;
                    var = signalData.variation;
                    phi = signalData.phi;

                }

            } catch (Exception ex) {
                //
            }
        }


        SignalData signalData = new SignalData();
        signalData.corr = maxCorr;
        signalData.variation = var;
        signalData.zScore = maxZ;
        signalData.phi = phi;
        return signalData;

    }

    public void show(double[] v) {

        for (int i = 0; i < v.length; i++)
            System.out.print(NumberFormatUtils.format(v[i]) + ", ");
        System.out.println();
    }

    public void show(int[] v) {

        for (int i = 0; i < v.length; i++)
            System.out.print(v[i] + ", ");
        System.out.println();
    }

    public void run() throws Exception {

        List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");

        for (String row : rows) {

            String pdbCode = row.substring(0, 4);
            String pdbChain = row.substring(5, 6);
            String strClass = row.split("\t")[1];
            if (strClass.equals("1"))
                strClass = "1";
            else
                strClass = "0";

            this.analyse(pdbCode, pdbChain, strClass);
        }

//        // read from francois
//        ProteinCSVReader csvReader = new ProteinCSVReader();
//        List<Row> rowsFrancois = csvReader.getData("data/tapo/FrancoisData1802.csv");
//        for (Row row : rowsFrancois) {
//
////            if (row.getByEyeTot() == 1) {
////                this.analyse(row.getPdbCode(), row.getPdbChain(), "1");
////
////            } else
//            if (row.getByEyeTot() == 0) {
//                this.analyse(row.getPdbCode(), row.getPdbChain(), "0");
//
//            }
//
//        }
//
//        List<String> rows = DataIO.readLines("data/evaluation/solenoid.data.TRs");
//
//        for (String row : rows) {
//            String pdbCode = row.substring(0, 4);
//            String pdbChain = row.substring(5, 6);
//            String strClass = "1";
//            String entry = pdbCode + "_" + pdbChain;
//
//            this.analyse(pdbCode, pdbChain, "1");
//
//        }

//        // read from RepeatsDB.
//        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");
//
//
//        //convert to rows
//        for (
//                RowRepeatDB r
//                : rowsRB)
//
//        {
//
//            try {
//                if (r.getAnnlevel().equals("Detailed") && !r.getStrclass().equals("II.2")) {
//                    int noTRs = r.getUnits().split(";").length;
//                    int startRegion = 0;
//                    int endRegion = 0;
//                    String region = r.getRegion();
//                    if (region.startsWith("-")) {
//                        region = region.substring(1, region.length());
//                        startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
//                        endRegion = Integer.parseInt(region.split("-")[1]);
//
//                    } else {
//                        startRegion = Integer.parseInt(region.split("-")[0]);
//                        endRegion = Integer.parseInt(region.split("-")[1]);
//                    }
//                    double avgL = (double) (endRegion - startRegion + 1) / noTRs;
//
//                    //if (avgL < 25) {
//
//                    String entry = r.getPdbCode() + "_" + r.getPdbChain();
//
//                    // do something here.
//
//                    this.analyse(r.getPdbCode(), r.getPdbChain(), "1");
//
//
//                    // }
//                }
//
//            } catch (Exception ex) {
//
//            }


//        }

    }


}
