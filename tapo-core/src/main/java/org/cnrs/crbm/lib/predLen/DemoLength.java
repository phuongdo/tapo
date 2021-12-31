package org.cnrs.crbm.lib.predLen;

import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.raphael.Raphael;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.maths.LocationWrapper;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 8/14/2015.
 */
public class DemoLength {

    ReadFasta readFasta = new ReadFasta();
    Map<String, TaPoFastaFormat> fastaTRs;

    public DemoLength(String dir) throws Exception {
        fastaTRs = readFasta.readTaPoFastaFormat(dir);
//        System.out.println(fastaTRs.size());
//        System.exit(0);
    }

    public DemoLength() {
    }

    public static void main(String[] args) {

        //
        try {
            new DemoLength("data/predLen/tapov1.1.0").test5PredictedByJustTAPO_OUTPUT("1a9n_C", 0, "");
//            new DemoLength("").benchmark();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void benchmark() throws Exception {
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        for (String row : releventTRs) {
            try {
                String[] data = row.split("\t");
                String pdb = data[0];
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);

                int startRegion = 0;
                int endRegion = 0;
                String region = data[2];
                String strClass = data[1];
                double avgRefL = Double.parseDouble(data[3]);

                this.test5PredictedByJustTAPO(pdb, avgRefL, region);
            } catch (Exception ex) {
//                ex.printStackTrace();
            }
        }
    }

    public void test5PredictedByJustTAPO_OUTPUT(String pdb, double refLen, String refRegion) throws Exception {
        //String pdb = "1fbl_A";
        TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);
        List<Repeat> lstRep = taPoFastaFormat.getRepeats();
        String pdbCode = pdb.substring(0, 4);
        String pdbChain = pdb.substring(5, 6);
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        Atom[] atoms = repeatFinder.getAtoms();
        String strSS = repeatFinder.getStrSS();
        double[] consensus = new double[atoms.length];
        double[] consensusSize = new double[atoms.length];
        double[] consensusLen = new double[atoms.length];
        for (int i = 0; i < consensusLen.length; i++) {
            consensusLen[i] = 1000;
        }
        for (int i = 0; i < consensusSize.length; i++) {
            int size = 0;
            for (Repeat repeat : lstRep) {
                if (repeat.getStart() <= i && i <= repeat.getEnd()) {
                    size++;
                }
            }
            consensusSize[i] = size;
        }
        Map<Integer, List<Double>> map = new HashMap<Integer, List<Double>>();
        for (Repeat repeat : lstRep) {
            for (int i = repeat.getStart(); i < repeat.getEnd(); i++) {
                consensus[i] += repeat.getScore();
                //consensusLen[i] = Math.min(consensusLen[i], repeat.getAvgLength());
            }
            for (RepeatContent rc : repeat.getRepeats()) {
                //RepeatContent rc = repeat.getRepeats().get(0);
                String pattern = VectorShape.getSSPattern(strSS
                        .substring(rc.getStart(), rc.getEnd()));
                for (int i = rc.getStart(); i < rc.getEnd(); i++) {

                    if (pattern.length() >= 2)
                        consensusLen[i] = Math.min(consensusLen[i], rc.size());
                }

                if (pattern.length() >= 2) {

                    if (map.containsKey(pattern.length())) {
                        List<Double> lstD = map.get(pattern.length());
                        lstD.add(repeat.getAvgLength());
                        map.put(pattern.length(), lstD);

                    } else {
                        List<Double> lstD = new ArrayList<Double>();
                        lstD.add(repeat.getAvgLength());
                        map.put(pattern.length(), lstD);

                    }


                }
            }


        }
        //normalize
        for (int i = 0; i < consensus.length; i++) {
            if (consensusSize[i] > 0)
                consensus[i] = consensus[i] / consensusSize[i];

        }

//        for (int i = 0; i < consensusLen.length; i++) {
//            System.out.print(NumberFormatUtils.format(consensus[i]) + ",");
//        }




        List<Region> lstRegions = this.findRegion(consensus);
        System.out.println(lstRegions);
        for (Region rc : lstRegions) {
//            System.out.println(rc);

            double predLen = 0.0;

            String predRegion = PdbTools.getResSeq(atoms, rc.getStart()) + "-" + PdbTools.getResSeq(atoms, rc.getEnd());

            double minL = 1000;
            for (Repeat repeat : lstRep) {
                if (repeat.getFinderName().equals("RMSD") || repeat.getFinderName().equals("VECTORS")) {
                    minL = Math.min(minL, repeat.getAvgLength());
                }
            }

            predLen = minL;

            if (predLen > 0) {
                System.out.println(pdb + "\t" + refRegion + "\t" + refLen + "\t" + predRegion + "\t" + NumberFormatUtils.format(predLen));
//                System.out.println(lstPreds);
            }
        }


    }


    public void test5PredictedByJustTAPO(String pdb, double refLen, String refRegion) throws Exception {

        //String pdb = "1fbl_A";
        TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);
        List<Repeat> lstRep = taPoFastaFormat.getRepeats();
        String pdbCode = pdb.substring(0, 4);
        String pdbChain = pdb.substring(5, 6);
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        Atom[] atoms = repeatFinder.getAtoms();
        String strSS = repeatFinder.getStrSS();
        double[] consensus = new double[atoms.length];
        double[] consensusSize = new double[atoms.length];
        double[] consensusLen = new double[atoms.length];

        for (int i = 0; i < consensusLen.length; i++) {
            consensusLen[i] = 1000;

        }
        for (int i = 0; i < consensusSize.length; i++) {

            int size = 0;
            for (Repeat repeat : lstRep) {
                if (repeat.getStart() < i && i < repeat.getEnd()) {
                    size++;
                }
            }
            consensusSize[i] = size;

        }
        Map<Integer, List<Double>> map = new HashMap<Integer, List<Double>>();
        for (Repeat repeat : lstRep) {

            for (int i = repeat.getStart(); i < repeat.getEnd(); i++) {
                consensus[i] += repeat.getScore();
                //consensusLen[i] = Math.min(consensusLen[i], repeat.getAvgLength());
            }

            for (RepeatContent rc : repeat.getRepeats()) {
                //RepeatContent rc = repeat.getRepeats().get(0);
                String pattern = VectorShape.getSSPattern(strSS
                        .substring(rc.getStart(), rc.getEnd()));

                for (int i = rc.getStart(); i < rc.getEnd(); i++) {

                    if (pattern.length() >= 2)
                        consensusLen[i] = Math.min(consensusLen[i], rc.size());
                }

                if (pattern.length() >= 2) {

                    if (map.containsKey(pattern.length())) {
                        List<Double> lstD = map.get(pattern.length());
                        lstD.add(repeat.getAvgLength());
                        map.put(pattern.length(), lstD);

                    } else {
                        List<Double> lstD = new ArrayList<Double>();
                        lstD.add(repeat.getAvgLength());
                        map.put(pattern.length(), lstD);

                    }


                }
            }


        }
        //normalize
        for (int i = 0; i < consensus.length; i++) {
            if (consensusSize[i] > 0)
                consensus[i] = consensus[i] / consensusSize[i];

        }

//        for (int i = 0; i < consensusLen.length; i++) {
//            System.out.print(NumberFormatUtils.format(consensus[i]) + ",");
//        }
//
//        System.exit(0);


        List<Region> lstRegions = this.findRegion(consensus);
        for (Region rc : lstRegions) {
//            System.out.println(rc);


            Atom[] fraAtoms = Fragement.getFragementsofAtoms(atoms, rc.getStart(), rc.getEnd());
            double predLenCE = this.getPredictedLengthByCE(fraAtoms);

            // predict by consensus
            double predLenCONSENCUS = 0.0;
            int index = 0;
            for (int i = rc.getStart(); i <= rc.getEnd(); i++) {
                if (consensusLen[i] < 1000) {
                    predLenCONSENCUS += consensusLen[i];
                    index++;
                }
            }
            if (index > 0)
                predLenCONSENCUS = predLenCONSENCUS / rc.size();

            // predict by Raphael
            Raphael raphael = new Raphael(fraAtoms);
            double predRaphael = raphael.getRepeatLength();

            List<Double> lstPreds = new ArrayList<Double>();
            if (predLenCE > 20)
                lstPreds.add(predLenCE);
            if (predRaphael > 20) {
                lstPreds.add(predRaphael);
            }
            if (predLenCONSENCUS > 0) {
                lstPreds.add(predLenCONSENCUS);
            }


            double predLen = 0.0;// finalLen(lstPreds.toArray(new Double[lstPreds.size()]));
            if (predRaphael - 5 <= predLenCE & predLenCE <= predRaphael + 5) {
                predLen = (predLenCE + predRaphael) / 2;
            } else {
                if (predLenCE > 20)
                    predLen = predLenCE;
                else if (predLenCE < 20 && predRaphael > 20 && raphael.getTotalScore() > 5) {
                    predLen = predRaphael;
                } else {
                    predLen = predLenCONSENCUS;
                }

            }

            String predRegion = PdbTools.getResSeq(atoms, rc.getStart()) + "-" + PdbTools.getResSeq(atoms, rc.getEnd());

            if (predLen > 0) {
                System.out.println(pdb + "\t" + refRegion + "\t" + refLen + "\t" + predRegion + "\t" + NumberFormatUtils.format(predLen) + "\t" + NumberFormatUtils.format(predictLength(lstPreds)));
//                System.out.println(lstPreds);
            }
        }

//        System.out.println();
//        int minMotif = 10000;
//        List<Double> lstD = new ArrayList<Double>();
//        for (Map.Entry<Integer, List<Double>> entry : map.entrySet()) {
//            if (entry.getKey() < minMotif) {
//                minMotif = entry.getKey();
//                lstD = entry.getValue();
//            }
//
//        }
//
//        double predLen = 0.0;
//        for (Double a : lstD) {
//            predLen += a;
//
//        }
//
//        if (predLen > 0) {
//            System.out.println("pred:" + predLen / lstD.size());
//        }


    }

    public double finalLen(Double[] predLens) {
        double predLen = 0.0;


        double T = 5;

        int[] labels = new int[predLens.length];
        for (int i = 0; i < labels.length; i++)
            labels[i] = -1;

        int label = 0;
        double minLen = 1000000;
        for (int i = 0; i < predLens.length; i++) {
            minLen = Math.min(minLen, predLens[i]);
            if (labels[i] == -1) {
                double periodRf = predLens[i];
                labels[i] = label;

                for (int j = i + 1; j < predLens.length; j++) {
                    if (labels[j] == -1)
                        if (periodRf - T <= predLens[j]
                                && predLens[j] <= periodRf + T) {
                            labels[j] = label;
                        }
                }

                label++;
            }

        }

        if (label == 3) {
            int size = 0;
            double avgLen = 0.0;
            for (int i = 0; i < labels.length; i++) {
                avgLen += predLens[i];
                size++;

            }
            if (size > 0)
                predLen = avgLen / size;
        } else if (label == 2) {

            Map<Integer, Integer> C = new Raphael().getCount(labels);

            int maxLab = 0;
            int indexLab = -1;
            for (int i = 0; i < label; i++) {

                int freq = C.get(i);
                if (freq > maxLab) {
                    maxLab = freq;
                    indexLab = i;
                }
            }

            int size = 0;
            double avgLen = 0.0;
            for (int i = 0; i < labels.length; i++) {
                if (labels[i] == indexLab) {
                    avgLen += predLens[i];
                    size++;
                }
            }
            if (size > 0)
                predLen = avgLen / size;


        } else if (label == 1) {
            double avgLen = 0.0;
            for (int i = 0; i < predLens.length; i++) {
                avgLen += predLens[i];
            }
            predLen = avgLen / predLens.length;

        }

        return predLen;
    }

    final double THRESHOLD_CUTOFF = 0.3;

    public List<Region> findRegion(double[] concescus) {
        List<Region> lstRegions = new ArrayList<Region>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < concescus.length; i++) {
            double ch = concescus[i];
            if (ch >= THRESHOLD_CUTOFF) {
                start = i;
                while (i < concescus.length
                        && concescus[i] >= THRESHOLD_CUTOFF) {
                    i++;
                }
                i = end = i - 1;
                // save
                RepeatContent rc = new RepeatContent(start, end);
                if (rc.size() > 10)
                    lstRegions.add(rc);
            }


        }
        return lstRegions;
    }

    public void test4PredictedByCEBased() throws Exception {
        //read the reference dataset
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/predLen/demoPredLen.txt");

        for (String row : releventTRs) {
            try {
                String[] data = row.split("\t");
                String pdb = data[0];
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);
//            System.out.println(pdb);
                int startRegion = 0;
                int endRegion = 0;
                String region = data[2];
                String strClass = data[1];
                double avgRefL = Double.parseDouble(data[3]);
//                if (region.startsWith("-")) {
//                    region = region.substring(1, region.length());
//                    startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
//                    endRegion = Integer.parseInt(region.split("-")[1]);
//
//                } else {
//                    startRegion = Integer.parseInt(region.split("-")[0]);
//                    endRegion = Integer.parseInt(region.split("-")[1]);
//                }
//
////                if (endRegion - startRegion > 600)
////                    continue;
//                RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
//                Atom[] atoms = repeatFinder.getAtoms();
//                String strSS = repeatFinder.getStrSS();
//                TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);
//                List<Repeat> lstRep = taPoFastaFormat.getRepeats();
//                // group repeat into cluster
//                Map<String, List<Repeat>> repeatsGroup = new HashMap<String, List<Repeat>>();
//                for (Repeat repeat : lstRep) {
//
//                    if (repeatsGroup.containsKey(repeat.getCluster())) {
//                        List<Repeat> group = repeatsGroup.get(repeat.getCluster());
//                        group.add(repeat);
//                        repeatsGroup.put(repeat.getCluster(), group);
//                    } else {
//                        List<Repeat> group = new ArrayList<Repeat>();
//                        group.add(repeat);
//                        repeatsGroup.put(repeat.getCluster(), group);
//                    }
//
//                }
//
////                if(pdb.equals("1ofl_A"))
////                    System.out.println();
//                // convert to cluster
//                for (Map.Entry<String, List<Repeat>> group : repeatsGroup.entrySet()) {
//                    // each op group;
//                    List<Repeat> lstRepGroup = group.getValue();
//                    List<Double> lstLens = new ArrayList<Double>();
//                    int regionStart = 1000000;
//                    int regionsEnd = -1;
//
//                    for (Repeat repeat : lstRepGroup) {
//
//                        regionStart = Math.min(repeat.getStart(), regionStart);
//                        regionsEnd = Math.max(repeat.getEnd(), regionsEnd);
//                    }

//                Atom[] fraAtoms = Fragement.getFragementsofAtoms(atoms, regionStart, regionsEnd);
//                predLen = this.getPredictedLengthByCE(fraAtoms);
//            }


                double predLen = 0.0;


//                if (predLen > 0)
//                    System.out.println(pdb + "\t" + group.getKey() + "\t" + startRegion + "-" + endRegion + "\t" + avgRefL + "\t" + regionStart + "-" + regionsEnd + "\t" + NumberFormatUtils.format(predLen));


            } catch (Exception ex) {
                //do something
                ex.printStackTrace();
            }
        }


    }

    public void test3() throws Exception {
        //read the reference dataset
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/predLen/demoPredLen.txt");

        for (String row : releventTRs) {


            try {
                String[] data = row.split("\t");
                String pdb = data[0];
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);
//            System.out.println(pdb);
                int startRegion = 0;
                int endRegion = 0;
                String region = data[2];
                String strClass = data[1];
                double avgRefL = Double.parseDouble(data[3]);
                if (region.startsWith("-")) {
                    region = region.substring(1, region.length());
                    startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                    endRegion = Integer.parseInt(region.split("-")[1]);

                } else {
                    startRegion = Integer.parseInt(region.split("-")[0]);
                    endRegion = Integer.parseInt(region.split("-")[1]);
                }

//                if (endRegion - startRegion > 600)
//                    continue;
                RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                Atom[] atoms = repeatFinder.getAtoms();
                String strSS = repeatFinder.getStrSS();


                TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);
                List<Repeat> lstRep = taPoFastaFormat.getRepeats();
                // group repeat into cluster
                Map<String, List<Repeat>> repeatsGroup = new HashMap<String, List<Repeat>>();
                for (Repeat repeat : lstRep) {

                    if (repeatsGroup.containsKey(repeat.getCluster())) {
                        List<Repeat> group = repeatsGroup.get(repeat.getCluster());
                        group.add(repeat);
                        repeatsGroup.put(repeat.getCluster(), group);
                    } else {
                        List<Repeat> group = new ArrayList<Repeat>();
                        group.add(repeat);
                        repeatsGroup.put(repeat.getCluster(), group);
                    }

                }

//                if(pdb.equals("1ofl_A"))
//                    System.out.println();
                // convert to cluster
                for (Map.Entry<String, List<Repeat>> group : repeatsGroup.entrySet()) {
                    // each op group;
                    List<Repeat> lstRepGroup = group.getValue();
                    List<Double> lstLens = new ArrayList<Double>();
                    int regionStart = 1000000;
                    int regionsEnd = -1;

                    for (Repeat repeat : lstRepGroup) {

                        double avgPredLen = repeat.getAvgLength();
//            System.out.println(repeat);
                        double tmpRredLen = 0;
                        for (RepeatContent rc : repeat.getRepeats()) {
                            // RepeatContent rc = repeat.getRepeats().get(0);
                            String pattern = VectorShape.getSSPattern(strSS
                                    .substring(rc.getStart(), rc.getEnd()));
//                System.out.println(pattern);
                            if (pattern.length() >= 3) {
//                    System.out.println(">>parse ");
                                Atom[] fraAtoms = Fragement.getFragementsofAtoms(repeatFinder.getAtoms(), rc.getStart(), rc.getEnd());
                                double avgLenByCE = this.getPredictedLengthByCE(fraAtoms);
                                if (avgLenByCE > 10) {
                                    //avgPredLen = Math.min(avgPredLen, avgLenByCE);
                                    lstLens.add(avgLenByCE);
                                    tmpRredLen += avgLenByCE;

                                } else
                                    tmpRredLen += rc.size();

//                                else
//                                    lstLens.add(avgPredLen);

//                            } else if (pattern.length() >= 2)
//                                lstLens.add(avgPredLen);
                            } else if (pattern.length() >= 2)
                                tmpRredLen += rc.size();


                        }

                        tmpRredLen = tmpRredLen / repeat.getRepeats().size();
                        lstLens.add(tmpRredLen);
                        regionStart = Math.min(repeat.getStart(), regionStart);
                        regionsEnd = Math.max(repeat.getEnd(), regionsEnd);
                    }
//                    System.out.println(pdb + "\t" + regionStart + "-" + regionsEnd + "\t" + NumberFormatUtils.format(avgLenByCE));
                    double predLen = 0.0;
                    if (lstLens.size() > 0) {
                        predLen = this.predictLength(lstLens);
                    }
//                    else {
//                        Atom[] fraAtoms = Fragement.getFragementsofAtoms(repeatFinder.getAtoms(), regionStart, regionsEnd);
//                        predLen = this.getPredictedLengthByCE(fraAtoms);
//
//                    }

                    if (predLen > 0)
                        System.out.println(pdb + "\t" + group.getKey() + "\t" + startRegion + "-" + endRegion + "\t" + avgRefL + "\t" + regionStart + "-" + regionsEnd + "\t" + NumberFormatUtils.format(predLen));

                }


            } catch (Exception ex) {
                //do something
                ex.printStackTrace();
            }
        }


    }

    public void test2() throws Exception {
        // test 2

        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/predLen/demoPredLen.txt");


        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {

            //String pdb = "1io0_A";
            String pdb = entry.getKey();

            String pdbCode = pdb.substring(0, 4);
            String pdbChain = pdb.substring(5, 6);
            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            Atom[] atoms = repeatFinder.getAtoms();
            String strSS = repeatFinder.getStrSS();

            TaPoFastaFormat taPoFastaFormat = entry.getValue(); //fastaTRs.get(pdb);
            List<Repeat> lstRep = taPoFastaFormat.getRepeats();
            // group repeat into cluster
            Map<String, List<Repeat>> repeatsGroup = new HashMap<String, List<Repeat>>();
            for (Repeat repeat : lstRep) {
                if (repeatsGroup.containsKey(repeat.getCluster())) {
                    List<Repeat> group = repeatsGroup.get(repeat.getCluster());
                    group.add(repeat);
                    repeatsGroup.put(repeat.getCluster(), group);
                } else {
                    List<Repeat> group = new ArrayList<Repeat>();
                    group.add(repeat);
                    repeatsGroup.put(repeat.getCluster(), group);
                }

            }


            // convert to cluster
            for (Map.Entry<String, List<Repeat>> group : repeatsGroup.entrySet()) {

                // each op group;

                List<Repeat> lstRepGroup = group.getValue();
                List<Double> lstLens = new ArrayList<Double>();
                int regionStart = 1000000;
                int regionsEnd = -1;

                for (Repeat repeat : lstRepGroup) {

//                    double avgPredLen = repeat.getAvgLength();
////            System.out.println(repeat);
//                    double tmpRredLen = 0;
//                    for (RepeatContent rc : repeat.getRepeats()) {
////            RepeatContent rc = repeat.getRepeats().get(0);
//                        String pattern = VectorShape.getSSPattern(strSS
//                                .substring(rc.getStart(), rc.getEnd()));
////                System.out.println(pattern);
//                        if (pattern.length() >= 4) {
////                    System.out.println(">>parse ");
//                            Atom[] fraAtoms = Fragement.getFragementsofAtoms(repeatFinder.getAtoms(), rc.getStart(), rc.getEnd());
//                            double avgLenByCE = this.getPredictedLengthByCE(fraAtoms);
//                            if (avgLenByCE > 10) {
//                                avgPredLen = Math.min(avgPredLen, avgLenByCE);
//                            }
//                        }
//                    }
//                    if (tmpRredLen > 10)
//                        avgPredLen = tmpRredLen;
////            System.out.println(avgPredLen);
//                    lstLens.add(avgPredLen);


                    regionStart = Math.min(repeat.getStart(), regionStart);
                    regionsEnd = Math.max(repeat.getEnd(), regionsEnd);
                }
                Atom[] fraAtoms = Fragement.getFragementsofAtoms(repeatFinder.getAtoms(), regionStart, regionsEnd);
                double avgLenByCE = this.getPredictedLengthByCE(fraAtoms);
                System.out.println(pdb + "\t" + regionStart + "-" + regionsEnd + "\t" + NumberFormatUtils.format(avgLenByCE));

            }


        }

    }


    public double predictLength(List<Double> lstLens) {


//
//        DescriptiveStatistics stats = new DescriptiveStatistics();
//        for (double avgL : lstLens) {
//            stats.addValue(avgL);
//        }
//        double P_avg = stats.getMean();
//        double P_sd = stats.getStandardDeviation();
//        stats = new DescriptiveStatistics();
//        for (double avgL : lstLens) {
//
//            if (P_avg - P_sd / 2 <= avgL
//                    && avgL <= P_avg + P_sd / 2)
//                stats.addValue(avgL);
//        }
//
//        return stats.getMean();

        double predAvgLen = 1000000.0;
        List<LocationWrapper> clusterInput = new ArrayList<LocationWrapper>();

        if (lstLens.size() >= 3) {
            for (double avgL : lstLens) {
                clusterInput.add(new LocationWrapper(avgL));
            }
            KMeansPlusPlusClusterer<LocationWrapper> clusterer = new KMeansPlusPlusClusterer<LocationWrapper>(
                    2, 10000);
            List<CentroidCluster<LocationWrapper>> clusterResults = clusterer
                    .cluster(clusterInput);


            double maxSize = 0.0;
            // output the clusters
            for (int i = 0; i < clusterResults.size(); i++) {
//            System.out.println(clusterResults.get(i).getCenter() + ":" + clusterResults.get(i).getPoints().size());
                if (maxSize < clusterResults.get(i).getPoints().size()) {
                    maxSize = clusterResults.get(i).getPoints().size();
                    predAvgLen = clusterResults.get(i).getCenter().getPoint()[0];
                }

            }
        } else {
            //smaller will be selected
            for (double avgL : lstLens) {
                predAvgLen = Math.min(avgL, predAvgLen);
            }
        }


        return predAvgLen;
    }

    public void test() throws Exception {

        //read the reference dataset
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        for (String row : releventTRs) {


            try {
                String[] data = row.split("\t");
                String pdb = data[0];
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);
//            System.out.println(pdb);
                int startRegion = 0;
                int endRegion = 0;
                String region = data[2];
                String strClass = data[1];
                double avgRefL = Double.parseDouble(data[3]);
                if (region.startsWith("-")) {
                    region = region.substring(1, region.length());
                    startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                    endRegion = Integer.parseInt(region.split("-")[1]);

                } else {
                    startRegion = Integer.parseInt(region.split("-")[0]);
                    endRegion = Integer.parseInt(region.split("-")[1]);
                }

                if (endRegion - startRegion > 600)
                    continue;

                RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                startRegion = PdbTools.getPosition(repeatFinder.getAtoms(), startRegion);
                endRegion = PdbTools.getPosition(repeatFinder.getAtoms(), endRegion);
                Atom[] atoms = Fragement.getFragementsofAtoms(repeatFinder.getAtoms(), startRegion, endRegion);
                Raphael raphael = new Raphael(atoms);
                double avgLenByCE = this.getPredictedLengthByCE(atoms);
                System.out.println(pdb + "\t" + avgRefL + "\t" + NumberFormatUtils.format(raphael.getRepeatLength()) + "\t" + NumberFormatUtils.format(avgLenByCE));

            } catch (Exception ex) {
                //do something
                ex.printStackTrace();
            }
        }

    }

    public double getPredictedLengthByCE(Atom[] atoms) throws Exception {

//        Atom[] ca1 = atoms;
//        Atom[] ca2 = StructureTools.cloneCAArray(ca1);
//        CeSymm ceSymm = new CeSymm();
//        AFPChain afpChain = ceSymm.pairAlign(ca1, ca2);
//        int symmNr = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
//        int[][][] optAln = afpChain.getOptAln();
//        int[] blockLen = afpChain.getOptLen();
//        // find the maximum block which is max length;
//        int indexBlock = 0;
//        int maxBlockLen = 0;
//        for (int block = 0; block < afpChain.getBlockNum(); block++) {
//            if (maxBlockLen < blockLen[block]) {
//                indexBlock = block;
//                maxBlockLen = blockLen[block];
//            }
//        }
//
//        double avgL = 0.0;
//        int noAlgRes = 0;
//        for (int i = 0; i < blockLen[indexBlock]; i++) {
//            int posA1 = optAln[indexBlock][0][i];
//            int posA2 = optAln[indexBlock][1][i];
//            avgL += Math.abs(posA1 - posA2);
//            noAlgRes++;
//
//        }
//
//
//        double L = 0;
//        if (symmNr >= 2)
//            L = atoms.length / symmNr;
//        else {
//            if (noAlgRes > 0 && (afpChain.getTMScore() > 0.4))
//                L = avgL / noAlgRes;
//        }
//
//        if ((symmNr + afpChain.getTMScore() < 1.4) || (L > atoms.length / 2) || (L > noAlgRes / 2))
//            return 0;
//
//
//        return L;
        return 0;
    }
}
