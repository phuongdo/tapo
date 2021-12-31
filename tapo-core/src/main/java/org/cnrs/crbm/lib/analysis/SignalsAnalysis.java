package org.cnrs.crbm.lib.analysis;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.math.*;
import org.cnrs.crbm.lib.Rlib.RDir;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import rcaller.RCaller;

import java.util.*;

/**
 * Created by pdoviet on 10/3/2014.
 */
public class SignalsAnalysis {


    public static void main(String[] args) throws Exception {

//         new SignalsAnalysis().run();


//        int[] his = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
//
//
//        int winsize = 6;
//
//        for (int i = 0; i < his.length; ) {
//            int start = i;
//            int end = i + winsize;
//
//            if (end > his.length)
//                end = his.length;
//            int[] v = Arrays.copyOfRange(his, start, end);
//
//            for (int j = 0; j < v.length; j++) {
//
//                System.out.print(v[j] + ";");
//            }
//            System.out.println();
//
//            i = i + winsize;
//
//        }


        //new SignalsAnalysis().findPeriodicity("1c1z", "A");
        new SignalsAnalysis().signalAnalysis();
    }




    public void repeatsDBAnalysis() throws Exception {


        ProteinCSVReader csvReader = new ProteinCSVReader();

        List<RowRepeatDB> rows = csvReader.getRepeatDB("data/RDB-dataset.tab");

        for (RowRepeatDB row : rows) {

            //

            if (row.getAnnlevel().equals("Detailed") && row.getUnits().length() > 0) {
                //System.out.println(row.getEntry());


                String[] units = row.getUnits().split(";");
                double avgLeng = 0.0;
                for (String unit : units) {
                    try {
                        int starti = Integer.parseInt(unit.split("-")[0]);
                        int endi = Integer.parseInt(unit.split("-")[1]);
                        avgLeng += endi - starti + 1;
                    } catch (Exception ex) {

                    }
                }

                avgLeng = avgLeng / units.length;

                System.out.println(row.getEntry() + "\t" + row.getStrclass() + " \t" + NumberFormatUtils.format(avgLeng));

                // here is predicted


                this.predictPdb(row.getPdbCode(), row.getPdbChain(), avgLeng);
            }
        }
    }


    public List<SignalPeriod> predictPdb(String pdbCode, String pdbName, double avgLeng) throws Exception {
        SignalProcess signalProcess = new SignalProcess();

        List<SignalPeriod> posibleLengs = new ArrayList<SignalPeriod>();
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbName);
        double[] his = repeatFinder.getCmHistos();
        int[] wlist = new int[]{60, 90, 120, 160};
        for (Integer winsize : wlist) {

            for (int i = 0; i < his.length; ) {

                try {
                    int start = i;
                    int end = i + winsize;

                    if (end > his.length)
                        end = his.length;
                    double[] v = Arrays.copyOfRange(his, start, end);

                    SignalPeriod signalPeriod = new SignalPeriod(start, end, v);
                    // process v here
                    signalPeriod = signalProcess.costanceDetect(signalPeriod);
                    double period = -1;
                    if (signalPeriod.getPeriod() == -1) {
                        signalPeriod = signalProcess.periodDetect(signalPeriod);

                    }


                    if (signalPeriod.getPeriod() > 10) {

                        if (avgLeng - 2 < signalPeriod.getPeriod() && signalPeriod.getPeriod() < avgLeng + 2)
                            System.err.println("predicted >>> " + signalPeriod + " : " + NumberFormatUtils.format(avgLeng));
                        else
                            System.out.println("predicted >>> " + signalPeriod + " : " + NumberFormatUtils.format(avgLeng));
                    }

                    i = i + winsize;
                } catch (Exception ex) {
                    ///
                    i = i + winsize;
                }

            }
        }

        return posibleLengs;
    }


    public void findPeriodicity(String pdbCode, String pdbName) throws Exception {
//        String pdbCode = "1cvm";
//        String pdbName = "A";


        Superimposer superimposer = new Superimposer();
        //MutilAlign superimposer = new MutilAlign();
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbName);

        Atom[] atoms = repeatFinder.getAtoms();
        List signals = new ArrayList<Double>();
        int winsize = 20;


        List<Double> listPeriods = new ArrayList<Double>();
        for (int t = 0; t < atoms.length - winsize; t = t + winsize) {


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
            double[] hisFilter = Filter.filter(his);
            PeriodFinder pfinder = new PeriodFinder();
            //System.out.println(fft.findPeriodicity(hisFilter));

            //FFTperiodicity fft = new FFTperiodicity();
            //System.out.println(fft.findPeriodicity(hisFilter));

            //double finalPeriodicity = pfinder.findPeriodicity(hisFilter);
            //
            listPeriods.add(pfinder.findPeriodicity(hisFilter));
            //listPeriods.add(fft.findPeriodicity(hisFilter));


        }


        List<Double> posibleLengs = new ArrayList<Double>();
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

        for (Map.Entry<Integer, Double> entry : C.entrySet()) {
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
        for (double l : posibleLengs) {
            System.out.println(l);
        }


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

        for (Map.Entry<Integer, Double> entry : avg.entrySet()) {
            int label = entry.getKey();

            if (count.containsKey(label)) {
                if (count.get(label) > 3)
                    avg.put(label, (double) avg.get(label) / count.get(label));
            }

        }

        return avg;
    }


    public void signalAnalysis() throws Exception {
        String pdbCode = "1gvm";
        String pdbName = "A";


        StringBuffer code = new StringBuffer();

        Superimposer superimposer = new Superimposer();
        //MutilAlign superimposer = new MutilAlign();
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbName);

        Atom[] atoms = repeatFinder.getAtoms();
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

                // AFPChain afpChain = superimposer.pairAlign(seedAtoms, compareAtoms);
                // rmsd = afpChain.getAlignScore();
                //rmsd = afpChain.getChainRmsd();
                // double tmScore = afpChain.getTMScore();
                //  int a = afpChain.getOptLength();
                //System.out.println(a + ":" + tmScore);

                signals.add(rmsd);


            } catch (StructureException e) {
                e.printStackTrace();
            }

        }


        double[] his = ArrayUtils.toPrimitive((Double[]) signals.toArray(new Double[signals.size()]));
        //double[] his = repeatFinder.getCmHistos();
        FFTperiodicity fft = new FFTperiodicity();

        //his = Filter.filter(his);
        fft.show(his);
        System.out.println("");
        fft.show(Filter.filter(his));

        System.exit(0);
        //procesSignalRange(his, 30, 40);
        code.append("y<-c(");
        for (int i = 0; i < his.length; i++) {
            if (i < his.length - 1)
                code.append(his[i] + ",");
            else
                code.append(his[i]);
        }
        code.append(");");
        double[] hisFilter = Filter.filter(his);

        //double[] hisFilter = new double[his.length];

        //hisFilter = Filter.fftAutoCorrelation(his);
        StringBuffer buffer = new StringBuffer();


        //PeriodFinder pfinder = new PeriodFinder();
        //System.out.println(fft.findPeriodicity(hisFilter));

        //double finalPeriodicity = pfinder.findPeriodicity(hisFilter);
        double finalPeriodicity = fft.findPeriodicity(hisFilter);
        code.append("x<-c(");
        for (int i = 0; i < hisFilter.length; i++) {
            if (i < hisFilter.length - 1) {
                code.append(hisFilter[i] + ",");
                buffer.append(NumberFormatUtils.format(hisFilter[i]) + " ");

            } else {
                code.append(hisFilter[i]);
                buffer.append(NumberFormatUtils.format(hisFilter[i]));
            }
        }
        code.append(");");
        code.append("png('C:/Users/pdoviet/Desktop/file.png');");
        code.append("plot.ts(y, col=\"grey\", xlab=\"periodicity " + finalPeriodicity + "\", ylab=\"y\");");
        code.append("lines(x, col=\"blue\", xlab=\"pos\", ylab=\"x\");");

        callR(code);
        //System.out.println(buffer.toString());


    }

    public void run() throws Exception {

        String pdbCode = "1gvm";
        String pdbName = "A";


        StringBuffer code = new StringBuffer();


        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbName);


        double[] his = repeatFinder.getCmHistos();
        //procesSignalRange(his, 30, 40);
        code.append("y<-c(");
        for (int i = 0; i < his.length; i++) {
            if (i < his.length - 1)
                code.append(his[i] + ",");
            else
                code.append(his[i]);
        }
        code.append(");");
        double[] hisFilter = Filter.filter(repeatFinder.getCmHistos());
        code.append("x<-c(");
        for (int i = 0; i < hisFilter.length; i++) {
            if (i < hisFilter.length - 1)
                code.append(hisFilter[i] + ",");
            else
                code.append(hisFilter[i]);
        }
        code.append(");");
        code.append("png('C:/Users/pdoviet/Desktop/file.png');");
        code.append("plot.ts(y, col=\"grey\", xlab=\"pos\", ylab=\"maximun distance\");");
        code.append("lines(x, col=\"blue\", xlab=\"pos\", ylab=\"maximun distance\");");

        callR(code);


        int[] wlist = new int[]{60, 90, 120, 160};
        // split the signals


        for (Integer winsize : wlist) {

            for (int i = 0; i < his.length; ) {
                int start = i;
                int end = i + winsize;

                if (end > his.length)
                    end = his.length;
                double[] v = Arrays.copyOfRange(his, start, end);

                // process v here


                i = i + winsize;

            }
        }

//        FFTperiodicity fft = new FFTperiodicity();
//        List<Double> posibleLengs = fft.findPeriodicityWithWindow(hisFilter,
//                winsize);
//
//        //List<Double> posibleLengs = processSignals(his);
//        System.out.println("Fourier...");
//        for (double l : posibleLengs) {
//            System.out.println(l);
//        }
//        System.out.println("Costance...");
//        posibleLengs = processSignals(his);
//        for (double l : posibleLengs) {
//            System.out.println(l);
//        }
    }

    public List<Double> processSignals(double[] his) {

        List<Double> posibleLengs = new ArrayList<Double>();
        posibleLengs.add(procesSignalRange(his, 40, 50));
        posibleLengs.add(procesSignalRange(his, 30, 40));
        posibleLengs.add(procesSignalRange(his, 20, 30));
        posibleLengs.add(procesSignalRange(his, 10, 20));
        return posibleLengs;

    }

    public double procesSignalRange(double[] his, double min, double max) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = 0; i < his.length; i++) {
            if (min < his[i] && his[i] < max) {
                stats.addValue(his[i]);
            } else {
                // his[i] = 0;
            }
        }
        double P_avg = stats.getMean();
        double P_sd = stats.getStandardDeviation();
        stats = new DescriptiveStatistics();
        for (int i = 0; i < his.length; i++) {
            if (min < his[i] && his[i] < max) {
                if (P_avg - P_sd / 2 <= his[i]
                        && his[i] <= P_avg + P_sd / 2)
                    stats.addValue(his[i]);
                stats.addValue(his[i]);
            } else {
                // his[i] = 0;
            }
        }


        return stats.getMean();
    }

    public void callR(StringBuffer code) {
        RCaller caller = new RCaller();

        // Path of your R localisation
        caller.setRScriptExecutableFile(RDir.EXECUTABLE_FILE);

        try {
            caller.RunRCode(code.toString(), false, false);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }


}
