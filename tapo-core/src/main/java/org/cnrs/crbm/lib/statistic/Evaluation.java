package org.cnrs.crbm.lib.statistic;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.analysis.ScoreAnalysis;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.repeats.clusters.ClusterRepeat;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by pdoviet on 1/27/2015.
 */
public class Evaluation {


    static double QA_SCORE_THRESHOLD = 0.5;

    public static void main(String[] args) throws Exception {
        Evaluation evaluation = new Evaluation();
        //evaluation.evaluate();
        //evaluation.evaluate2();
//         evaluation.evaluateSonsole();
        // evaluation.evaluateTapoWithConsoleDataset();

        //evaluation.evaluateRepeatsDB();

        //evaluation.evaluateQAScore();

        //evaluation.evaluationTMScore();
        //evaluation.showDavros();
        // evaluation.evaluateTapoWithRepeatsDB();

        //evaluation.evaluateQAScore();


//        List<String> lines = DataIO.readLines("data/davros/top100pdbs.in");
//        for (String line : lines) {
//            System.out.println(line + "\t1");
//
//        }

        // read from RepeatsDB.


        //  evaluation.evaluateTAPODataset4();

        evaluation.evaluateTapoWithBenmarkSet();
        //evaluation.convertDataFromPDBNumberToArrayNumber();

//        evaluation.correlationBetweenPredictedAndReference();

//        evaluation.bechmarkPredictedRepeatUnitLength();

    }


    public void bechmarkPredictedRepeatUnitLength() throws Exception {

        double THRESHOLD = 0.2;
        Map<String, List<Repeat>> retrievedTRs = new ScoreAnalysis().readResult("data/benchmarkset/TRsRegion.o");
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        List<String> svmSayNo = DataIO.readLines("data/benchmarkset/SVMSayNo.o");
        int noRetrievedTRs = DataIO.readLines("data/benchmarkset/TRsRegion.o").size();
        int noReleventTRs = releventTRs.size();

        Set<String> positiveSet = new HashSet<String>();
        //  Map<String,String > mapReferenceTRs = new HashMap<String, String>();
        int noFoundTRs = 0;
        ClusterLocation clusterLocation = new ClusterLocation();
        System.out.println("====PRINT OUTPUT====");
        // check ability of find
        for (String row : releventTRs) {

            String[] data = row.split("\t");
            String pdb = data[0];

//            System.out.println(pdb);
            positiveSet.add(pdb);
//            if(pdb.equals("3m7h_A"))
//                System.out.println();
            int startRegion = 0;
            int endRegion = 0;
            String region = data[2];
            String strClass = data[1];
            double avgRefL = Double.parseDouble(data[3]);
            if (strClass.startsWith("II."))
                avgRefL = 7;

            if (region.startsWith("-")) {
                region = region.substring(1, region.length());
                startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                endRegion = Integer.parseInt(region.split("-")[1]);

            } else {
                startRegion = Integer.parseInt(region.split("-")[0]);
                endRegion = Integer.parseInt(region.split("-")[1]);
            }

            // check found or not{
            // }
            boolean found = false;
            int noFound = 0;
            int noCorrect = 0;
            if (retrievedTRs.containsKey(pdb)) {
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);
                RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                Atom[] atoms = repeatFinder.getAtoms();
//                System.out.println(atoms.length);

                List<Repeat> allTRs = retrievedTRs.get(pdb);
                noFound = allTRs.size();

//                Repeat top = allTRs.get(0);
                for (Repeat top : allTRs) {

                    int endTop = top.getEnd();
                    if (endTop >= atoms.length)
                        endTop = atoms.length - 1;
                    int startTop = top.getStart();
                    if (startTop >= atoms.length)
                        startTop = atoms.length - 1;
                    int relativeStart = this.getResSeq(startTop, atoms);
                    int relativeEnd = this.getResSeq(endTop, atoms);
                    double avgLegPredicted = top.getAvgLength();

                    double overlap = Math.abs(avgLegPredicted - avgRefL) / avgRefL;

                    if (overlap < THRESHOLD) {
                        found = true;
                        //break;
                        noCorrect++;
                    } else {

                    }
//                System.out.println(strClass + "\t" + avgRefL + "\t" + avgLegPredicted);
                    break;

                }
            }
            if (found) {
                noFoundTRs++;
                //check svm
                if (svmSayNo.contains(pdb))
                    System.out.println(row + "\t" + 0 + "\t0\t0");
                else
                    System.out.println(row + "\t" + 1 + "\t" + noCorrect + "\t" + noFound);

            } else {
                System.out.println(row + "\t" + 0 + "\t0\t0");
            }

        }


    }

    public void correlationBetweenPredictedAndReference() throws Exception {

        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        //

        Set<String> repeatDBs = new HashSet<String>();

        //convert to rows
        for (
                RowRepeatDB r
                : rowsRB)

        {

            try {
                if (r.getAnnlevel().equals("Detailed")) {
                    String entry = r.getPdbCode() + "_" + r.getPdbChain();
                    String region = r.getRegion();
                    repeatDBs.add(entry);
                    try {
                        if (r.getUnits().length() > 0) {
                            int noTRs = r.getUnits().split(";").length;
                            int startRegion = 0;
                            int endRegion = 0;

                            if (region.startsWith("-")) {
                                region = region.substring(1, region.length());
                                startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                                endRegion = Integer.parseInt(region.split("-")[1]);

                            } else {
                                startRegion = Integer.parseInt(region.split("-")[0]);
                                endRegion = Integer.parseInt(region.split("-")[1]);
                            }
                            double avgL = (double) (endRegion - startRegion + 1) / noTRs;

                            //if (avgL < 25) {


                            // do something here.
                            System.out.println(entry + "\t" + r.getStrclass() + "\t" + region + "\t" + NumberFormatUtils.format(avgL) + "\t" + noTRs);
                        } else {
                            System.out.println(entry + "\t" + r.getStrclass() + "\t" + region + "\t" + "1" + "\t" + 1);
                        }

                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }


                    // }
                }

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }


        Set<String> setFrancois = new HashSet<String>();
        List<String> francois = DataIO.readLines("data/benchmarkset/Francois.tmp");
        for (String row : francois) {
            String[] data = row.split("\t");
            String pdb = data[0];
            setFrancois.add(pdb);
            if (!repeatDBs.contains(pdb)) {
                System.out.println(pdb + "\t" + data[5] + "\t" + data[3] + "-" + data[4] + "\t" + data[1] + "\t" + data[2]);
            }

        }


        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        for (String row : releventTRs) {
            String[] data = row.split("\t");
            String pdb = data[0];
            if (!repeatDBs.contains(pdb) && !setFrancois.contains(pdb)) {
                System.out.println(row);
            }

        }

    }

    public String convertToRepeatDBClass(String strClass) {
        String out = "UNK";
        if (strClass.equals("21")) {
            out = "II.1";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        } else if (strClass.equals("22")) {
            out = "II.2";
        }


        return out;
    }

    public void convertDataFromPDBNumberToArrayNumber() {

        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.in");
        for (String row : releventTRs) {
            String[] data = row.split("\t");
            String pdb = data[0];


            int startRegion = 0;
            int endRegion = 0;
            String region = data[2];
            String strClass = data[1];
            if (region.startsWith("-")) {
                region = region.substring(1, region.length());
                startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                endRegion = Integer.parseInt(region.split("-")[1]);

            } else {
                startRegion = Integer.parseInt(region.split("-")[0]);
                endRegion = Integer.parseInt(region.split("-")[1]);
            }

            RepeatFinder repeatFinder = new RepeatFinder(pdb.substring(0, 4), pdb.substring(5, 6));
            Atom[] atoms = repeatFinder.getAtoms();

            System.out.println(pdb + "\t" + strClass + "\t" + this.getPosition(atoms, startRegion) + "-" + this.getPosition(atoms, endRegion));

        }

    }

    int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {
            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;
        }
        return 0;
    }

    public void evaluateTapoWithBenmarkSet() throws Exception {

        double THRESHOLD = 0.2;
        Map<String, List<Repeat>> retrievedTRs = new ScoreAnalysis().readResult("data/benchmarkset/TRsRegion.o");
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        List<String> svmSayNo = DataIO.readLines("data/benchmarkset/SVMSayNo.o");
        int noRetrievedTRs = DataIO.readLines("data/benchmarkset/TRsRegion.o").size();
        int noReleventTRs = releventTRs.size();

        Set<String> positiveSet = new HashSet<String>();
        //  Map<String,String > mapReferenceTRs = new HashMap<String, String>();
        int noFoundTRs = 0;
        ClusterLocation clusterLocation = new ClusterLocation();
        System.out.println("====PRINT OUTPUT====");
        // check ability of find
        for (String row : releventTRs) {

            String[] data = row.split("\t");
            String pdb = data[0];

//            System.out.println(pdb);
            positiveSet.add(pdb);
//            if(pdb.equals("3m7h_A"))
//                System.out.println();
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

            // check found or not{
            // }
            boolean found = false;
            int noFound = 0;
            int noCorrect = 0;
            if (retrievedTRs.containsKey(pdb)) {
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);
                RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                Atom[] atoms = repeatFinder.getAtoms();
//                System.out.println(atoms.length);

                List<Repeat> allTRs = retrievedTRs.get(pdb);
                noFound = allTRs.size();

                //Repeat top = allTRs.get(0);
                for (Repeat top : allTRs) {

                    int endTop = top.getEnd();
                    if (endTop >= atoms.length)
                        endTop = atoms.length - 1;
                    int startTop = top.getStart();
                    if (startTop >= atoms.length)
                        startTop = atoms.length - 1;
                    int relativeStart = this.getResSeq(startTop, atoms);
                    int relativeEnd = this.getResSeq(endTop, atoms);
                    double avgLegPredicted = top.getAvgLength();


                    Repeat trsFound = new Repeat();
                    trsFound.getRepeats().add(new RepeatContent(relativeStart, relativeEnd));
                    Repeat trsReferences = new Repeat();
                    trsReferences.getRepeats().add(new RepeatContent(startRegion, endRegion));

//
//                    double overlap = clusterLocation.cover(trsFound, trsReferences);
//
//                    if (overlap > THRESHOLD) {
//                        found = true;
//                        //break;
//                        noCorrect++;
//                    } else {
//
//                    }
//                    System.out.println(strClass + "\t" + avgRefL + "\t" + avgLegPredicted);

                }
            }
            if (found) {
                noFoundTRs++;
                //check svm
                if (svmSayNo.contains(pdb))
                    System.out.println(row + "\t" + 0 + "\t0\t0");
                else
                    System.out.println(row + "\t" + 1 + "\t" + noCorrect + "\t" + noFound);

            } else {
                System.out.println(row + "\t" + 0 + "\t0\t0");
            }

        }

        //
        int noTruePositive = 0;
        int noFalsePositive = 0;
        System.out.println("===========CALCULATION OF PRECISION AND RECALL==============");
        for (Map.Entry<String, List<Repeat>> entry : retrievedTRs.entrySet()) {
            String key = entry.getKey();
            List<Repeat> allTRs = entry.getValue();

            if (!positiveSet.contains(key)) {
                noFalsePositive++;
//                System.out.println(key);
            }

            for (Repeat top : allTRs) {
                for (String row : releventTRs) {
                    String[] data = row.split("\t");
                    String pdb = data[0];
                    if (pdb.equals(key)) {

                        int startRegion = 0;
                        int endRegion = 0;
                        String region = data[2];
                        String strClass = data[1];
                        if (region.startsWith("-")) {
                            region = region.substring(1, region.length());
                            startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                            endRegion = Integer.parseInt(region.split("-")[1]);

                        } else {
                            startRegion = Integer.parseInt(region.split("-")[0]);
                            endRegion = Integer.parseInt(region.split("-")[1]);
                        }
                        // check found or not{
                        // }

                        Repeat repeatXX = new Repeat();
                        repeatXX.getRepeats().add(new RepeatContent(startRegion, endRegion));
//                        double overlap = clusterLocation.cover(top, repeatXX);
//                        if (overlap > THRESHOLD) {
//                            noTruePositive++;
//                        }
                    }


                }

            }

        }


        double precision = (double) noTruePositive / noRetrievedTRs;
        double recall = (double) noFoundTRs / noReleventTRs;
        double sensitivity = recall;
        double specificity = (double) (132 - noFalsePositive) / 132;
        System.out.println("precision: " + NumberFormatUtils.format(precision));
        System.out.println("recall: " + NumberFormatUtils.format(recall));
        System.out.println("sensitivity: " + NumberFormatUtils.format(sensitivity));
        System.out.println("specificity: " + NumberFormatUtils.format(specificity));


    }


    private int getResSeq(int pos, Atom[] atoms) {

        return atoms[pos].getGroup().getResidueNumber().getSeqNum();
    }


    public void evaluateTAPODataset4() throws Exception {


        Map<String, List<Repeat>> map = new ScoreAnalysis().readResult("data/tapo/dataset4.out");
        List<String> rows = DataIO.readLines("data/benchmarkset/dataset4.in");

        for (String row : rows) {
            if (row.startsWith("#"))
                continue;


            String[] datas = row.split("\t");

            if (map.containsKey(datas[0])) {

                Repeat repeat = map.get(datas[0]).get(0);
                System.out.println(row + "\t1\t" + repeat.getScore());
            } else
                System.out.println(row + "\t0\t0");

        }

    }


    public void evaluateTapoWithRepeatsDB() throws Exception {
        Map<String, List<Repeat>> map = new ScoreAnalysis().readResult("data/tapo/result_Dataset5.o");
        ClusterLocation clusterLocation = new ClusterLocation();


        // read from RepeatsDB.

        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");


        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            try {
                if (r.getAnnlevel().equals("Detailed")) {
                    int noTRs = r.getUnits().split(";").length;
                    int startRegion = 0;
                    int endRegion = 0;
                    String region = r.getRegion();
                    if (region.startsWith("-")) {
                        region = region.substring(1, region.length());
                        startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);

                    } else {
                        startRegion = Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);
                    }
                    double avgL = (double) (endRegion - startRegion + 1) / noTRs;

                    String proteinCode = r.getPdbCode() + "_" + r.getPdbChain();

                    int count = 0;
                    if (map.containsKey(proteinCode)) {
                        List<Repeat> allTRs = map.get(proteinCode);
                        List<ClusterRepeat> clusters = clusterLocation.cluster(allTRs);

                        boolean found = false;
                        for (ClusterRepeat c : clusters) {
                            StringBuffer contentBuilder = new StringBuffer();
                            // System.out.println("-----");
                            c.sortScoreDESC();
                            //get top
                            Repeat top = c.getRepeats().get(0);
                            int start = top.getStart();
                            int end = top.getEnd();
                            // pdbId_chain Finder avgLeng nUnits RL
                            StringBuffer unitsBuffer = new StringBuffer();
                            for (RepeatContent unit : top.getRepeats()) {
                                unitsBuffer.append(unit.getStart() + "-" + unit.getEnd() + ";");
                            }
                            String unitsStr = unitsBuffer.toString();
                            unitsStr = unitsStr.substring(0, unitsStr.length() - 1);
                            contentBuilder.append(
                                    +top.getRepeats().size() + "\t"
                                            + NumberFormatUtils.format(top.getAvgLength()) + "\t" + start
                                            + "\t" + end + "\t" + unitsStr + "\t" + NumberFormatUtils.format(top.getScore()) + "\t" + NumberFormatUtils.format(top.getRankScore()))
                            ;

                            Repeat repeatXX = new Repeat();
                            repeatXX.getRepeats().add(new RepeatContent(startRegion, endRegion));
//                            double overlap = clusterLocation.cover(top, repeatXX);
//
//
//                            if (overlap > 0.5) {
//                                System.out.println(proteinCode + "\t" + r.getStrclass() + "\t1" + "\t" + startRegion + "\t" + endRegion + "\t" + contentBuilder.toString() + "\t" + NumberFormatUtils.format(overlap) + "\t" + noTRs);
//                                found = true;
//                                count++;
//                            }

                        }
                        if (count > 1)
                            System.out.println(">>before");

                        if (!found)
                            System.out.println(proteinCode + "\t" + r.getStrclass() + "\t0");

                    } else {
                        System.out.println(proteinCode + "\t" + r.getStrclass() + "\t0");
                    }


                }

            } catch (Exception ex) {

            }


        }


        System.out.println("===========NO TR Class============");

        // read data
        csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");
        for (Row row : rowsFrancois) {
            if (row.getByEyeTot() == 2) {
                if (map.containsKey(row.getProtein())) {
                    System.out.println(row.getProtein() + "\t1");
                } else
                    System.out.println(row.getProtein() + "\t0");
            }

        }


    }

    public void evaluationTMScore() throws Exception {

        List<String> lines = DataIO.readLines("data/train/tmScore-10New.o");
        Map<String, String> map = new HashMap<String, String>();

        for (String line : lines) {
            //System.out.println(line);
            if (line.length() < 10)
                continue;
            String[] data = line.split("\t");
            map.put(data[0], data[4] + "\t" + data[1]);

        }

        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        List<Row> rows = new ArrayList<Row>();


        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            try {
                if (r.getAnnlevel().equals("Detailed")) {
                    int noTRs = r.getUnits().split(";").length;
                    int startRegion = 0;
                    int endRegion = 0;
                    String region = r.getRegion();
                    if (region.startsWith("-")) {
                        region = region.substring(1, region.length());
                        startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);

                    } else {
                        startRegion = Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);
                    }
                    double avgL = (double) (endRegion - startRegion + 1) / noTRs;

                    //if (avgL < 25) {

                    String entry = r.getPdbCode() + "_" + r.getPdbChain();

                    if (map.containsKey(entry)) {
                        System.out.println(entry + "\t" + map.get(entry) + "\t1");
                    }

                    // }
                }

            } catch (Exception ex) {

            }


        }

        csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");
        for (Row row : rowsFrancois) {
            if (row.getByEyeTot() == 2) {
                int classTRs = row.getByEyeTot();
                if (classTRs == 2)
                    classTRs = 0;
                String entry = row.getProtein();
                if (map.containsKey(entry)) {
                    System.out.println(entry + "\t" + map.get(entry) + "\t" + classTRs);
                }

            }

        }


    }


    public void evaluateTapoWith4Family() throws Exception {
        Map<String, List<Repeat>> map = new ScoreAnalysis().readResult("data/tapo/ExperimentalResult.o");


        List<String> lines = DataIO.readLines("data/davros/LRR.results");

        int numberMotifsFound = 0;
        int numberStructureFound = 0;
        for (String line : lines) {
            String[] data = line.split("\t");
            String protein = data[0];

            String proteinCode = "";
            if (protein.length() == 5)
                proteinCode = protein.substring(0, 4) + "_" + protein.substring(4, 5);
            else
                proteinCode = protein + "_A";


            if (map.containsKey(proteinCode)) {
                List<Repeat> allTRs = map.get(proteinCode);
                List<ClusterRepeat> clusters = new ClusterLocation().cluster(allTRs);
                ClusterLocation clusterLocation = new ClusterLocation();

                for (ClusterRepeat c : clusters) {
                    StringBuffer contentBuilder = new StringBuffer();
                    // System.out.println("-----");
                    c.sortScoreDESC();
                    //get top
                    Repeat top = c.getRepeats().get(0);
                    int start = top.getStart();
                    int end = top.getEnd();
                    // pdbId_chain Finder avgLeng nUnits RL
                    StringBuffer unitsBuffer = new StringBuffer();
                    for (RepeatContent unit : top.getRepeats()) {
                        unitsBuffer.append(unit.getStart() + "-" + unit.getEnd() + ";");
                    }
                    String unitsStr = unitsBuffer.toString();
                    unitsStr = unitsStr.substring(0, unitsStr.length() - 1);
                    contentBuilder.append(
                            +top.getRepeats().size() + "\t"
                                    + NumberFormatUtils.format(top.getAvgLength()) + "\t" + start
                                    + "\t" + end + "\t" + unitsStr + "\t" + NumberFormatUtils.format(top.getScore()) + "\t" + NumberFormatUtils.format(top.getRankScore()))
                    ;


                    System.out.println(proteinCode + "\t" + top.getStart() + "\t" + top.getEnd() + "\t" + contentBuilder.toString());


                }
            }


        }


    }


    public void showDavros() {

//        List<String> lines = DataIO.readLines("data/davros/top100_byTAPO.dat");
//
//        Map<String, String> map = new HashMap<String, String>();
//        for (String line : lines) {
//            String[] data = line.split("\t");
//            map.put(data[0], line);
//        }
//
//        for (Map.Entry<String, String> entry : map.entrySet()) {
//            System.out.println(entry.getKey());
//        }
//
//
        List<String> lines = DataIO.readLines("data/davros/aabarrel.results");
        for (String line : lines) {
            String[] data = line.split("\t");
            String protein = data[0];

            if (protein.length() == 5)
                System.out.println(protein.substring(0, 4) + "_" + protein.substring(4, 5));
            else
                System.out.println(protein + "_A");

        }
    }

    public void convert() {
        List<String> lines = DataIO.readLines("data/evaluation/solenoid.data.TRs");

        Map<String, String> map = new HashMap<String, String>();
        for (String line : lines) {
            String[] data = line.split("_");
            System.out.println(data[0].toLowerCase() + "_" + data[1]);

        }
    }

    public void evaluateTapoWithConsoleDataset() throws Exception {


        List<String> lines = DataIO.readLines("data/evaluation/tapo_console_set.o");
        Map<String, List<Data>> map = new HashMap<String, List<Data>>();
        for (String line : lines) {
            //System.out.println(line);
            if (line.length() < 10)
                continue;


            String[] data = line.split("\t");
//            System.out.println(line);
            int start = 0;
            int end = 0;
            try {
                if (data[5].startsWith("-")) {
                    data[5] = data[5].substring(1, data[5].length());
                    start = (-1) * Integer.parseInt(data[5].split("-")[0]);
                    end = Integer.parseInt(data[5].split("-")[1]);

                } else {
                    start = Integer.parseInt(data[5].split("-")[0]);
                    end = Integer.parseInt(data[5].split("-")[1]);
                }

            } catch (Exception ex) {
                ex.printStackTrace();
            }
            Data dataobj = new Data();
            dataobj.setPdbBegin(start);
            dataobj.setPdbEnd(end);
            dataobj.setPdbEntry(data[0]);
            dataobj.setScore1(Double.parseDouble(data[7]));
            dataobj.setScore2(Double.parseDouble(data[8]));

            if (map.containsKey(data[0])) {
                List<Data> l = map.get(dataobj.getPdbEntry());
                l.add(dataobj);
                map.put(dataobj.getPdbEntry(), l);
            } else {
                List<Data> l = new ArrayList<Data>();
                l.add(dataobj);
                map.put(dataobj.getPdbEntry(), l);
            }
        }


        // read data from console

        List<String> trs = DataIO.readLines("data/evaluation/solenoid.data.TRs");
        List<String> noTRs = DataIO.readLines("data/evaluation/solenoid.data.NoTRs");

        for (String tr : trs) {

            if (map.containsKey(tr)) {
                List<Data> l = map.get(tr);

                double maxQA = -1;
                double maxScore1 = -1;
                double maxScore2 = -1;
                for (Data data : l) {
                    maxQA = Math.max(maxQA, data.getQAScore());
                    maxScore1 = Math.max(maxScore1, data.getScore1());
                    maxScore2 = Math.max(maxScore2, data.getScore2());
                }
                if (maxQA >= QA_SCORE_THRESHOLD)
                    System.out.println(tr + "\t1\t1\t" + NumberFormatUtils.format(maxQA));
                else
                    System.out.println(tr + "\t1\t0\t" + NumberFormatUtils.format(maxQA));

            } else {
                System.out.println(tr + "\t1\t0\t0");
            }

        }

        for (String tr : noTRs) {
            if (map.containsKey(tr)) {

                List<Data> l = map.get(tr);

                double maxQA = -1;
                double maxScore1 = -1;
                double maxScore2 = -1;
                for (Data data : l) {
                    maxQA = Math.max(maxQA, data.getQAScore());
                    maxScore1 = Math.max(maxScore1, data.getScore1());
                    maxScore2 = Math.max(maxScore2, data.getScore2());
                }
                if (maxQA >= QA_SCORE_THRESHOLD)
                    System.out.println(tr + "\t0\t1\t" + NumberFormatUtils.format(maxQA));
                else
                    System.out.println(tr + "\t0\t0\t" + NumberFormatUtils.format(maxQA));


            } else {
                System.out.println(tr + "\t0\t0\t0");
            }

        }


    }

    public void evaluateSonsole() throws Exception {

        // read pdbList

        File folder = new File("output/console/converted");

        File[] listOfFiles = folder.listFiles();
        Set<String> set = new HashSet<String>();

        Map<String, String> map = new HashMap<String, String>();


        StringBuffer buffer = new StringBuffer();
        for (File file : listOfFiles) {
            if (file.isFile()) {
                //set.add(file.getName());
                //System.out.println(file);
                String[] data = DataIO.readFile(file.getAbsolutePath()).split("\t");
                map.put(data[0], data[4]);

            }
        }


        List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");

        for (String row : rows) {
            if (row.startsWith("#"))
                continue;

            String[] datas = row.split("\t");

            if (map.containsKey(datas[0])) {
                String a = map.get(datas[0]);
                a = a.substring(1, a.length() - 1);
                String b[] = a.split(" ");


                if (b.length > 0 && b[0].equals("True")) {
                    System.out.println(row + "\t1");
                } else
                    System.out.println(row + "\t0");

            } else
                System.out.println(row + "\t0");

        }


//        // read from francois
//        ProteinCSVReader csvReader = new ProteinCSVReader();
//        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");
//        for (Row row : rowsFrancois) {
//
//            if (row.getByEyeTot() == 1) {
//                if (map.containsKey(row.getProtein())) {
//
//                    if (map.get(row.getProtein()).contains("True"))
//                        System.out.println(row.getProtein() + "\t1\t1");
//                    else
//                        System.out.println(row.getProtein() + "\t1\t0");
//                }
//            } else if (row.getByEyeTot() == 2) {
//                if (map.containsKey(row.getProtein())) {
//                    if (map.get(row.getProtein()).contains("True"))
//                        System.out.println(row.getProtein() + "\t0\t1");
//                    else
//                        System.out.println(row.getProtein() + "\t0\t0");
//                }
//
//            }
//
//        }


    }

    public void evaluateQAScore() throws Exception {

        List<String> lines = DataIO.readLines("data/evaluation/QAEvaluation.o");

        Map<String, String> map = new HashMap<String, String>();
        for (String line : lines) {
            String[] data = line.split("\t");
            if (Double.parseDouble(data[1]) > 0)
                map.put(data[0], line);
        }


        // read from francois
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/tapo/FrancoisData1802.csv");
        for (Row row : rowsFrancois) {

            if (row.getByEyeTot() == 1) {

                if (map.containsKey(row.getProtein())) {
                    System.out.println(map.get(row.getProtein()) + "\tTRs");
                }
            } else if (row.getByEyeTot() == 0) {
                if (map.containsKey(row.getProtein())) {
                    System.out.println(map.get(row.getProtein()) + "\tNo-TRs");
                }

            }

        }

        // read from RepeatsDB.
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        List<Row> rows = new ArrayList<Row>();


        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            try {
                if (r.getAnnlevel().equals("Detailed")) {
                    int noTRs = r.getUnits().split(";").length;
                    int startRegion = 0;
                    int endRegion = 0;
                    String region = r.getRegion();
                    if (region.startsWith("-")) {
                        region = region.substring(1, region.length());
                        startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);

                    } else {
                        startRegion = Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);
                    }
                    double avgL = (double) (endRegion - startRegion + 1) / noTRs;

                    if (avgL < 25) {

                        String entry = r.getPdbCode() + "_" + r.getPdbChain();
                        if (map.containsKey(entry)) {
                            System.out.println(map.get(entry) + "\tTRs");
                        }
                    }
                }

            } catch (Exception ex) {

            }


        }


    }

    public void evaluateRepeatsDB() throws Exception {
        List<String> lines = DataIO.readLines("data/evaluation/RepeatsDB.o");
        Map<String, List<Data>> map = new HashMap<String, List<Data>>();
        for (String line : lines) {
            //System.out.println(line);
            if (line.length() < 10)
                continue;


            String[] data = line.split("\t");
            int start = 0;
            int end = 0;
            if (data[5].startsWith("-")) {
                data[5] = data[5].substring(1, data[5].length());
                start = (-1) * Integer.parseInt(data[5].split("-")[0]);
                end = Integer.parseInt(data[5].split("-")[1]);

            } else {
                start = Integer.parseInt(data[5].split("-")[0]);
                end = Integer.parseInt(data[5].split("-")[1]);
            }

            Data dataobj = new Data();
            dataobj.setPdbBegin(start);
            dataobj.setPdbEnd(end);
            dataobj.setPdbEntry(data[0]);


            if (map.containsKey(data[0])) {
                List<Data> l = map.get(dataobj.getPdbEntry());
                l.add(dataobj);
                map.put(dataobj.getPdbEntry(), l);
            } else {
                List<Data> l = new ArrayList<Data>();
                l.add(dataobj);
                map.put(dataobj.getPdbEntry(), l);
            }
        }

        System.out.println(map.size());


    }


    public void showTRsSmaller25Res() throws Exception {
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        List<Row> rows = new ArrayList<Row>();


        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            try {
                if (r.getAnnlevel().equals("Detailed")) {
                    int noTRs = r.getUnits().split(";").length;
                    int startRegion = 0;
                    int endRegion = 0;
                    String region = r.getRegion();
                    if (region.startsWith("-")) {
                        region = region.substring(1, region.length());
                        startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);

                    } else {
                        startRegion = Integer.parseInt(region.split("-")[0]);
                        endRegion = Integer.parseInt(region.split("-")[1]);
                    }
                    double avgL = (double) (endRegion - startRegion + 1) / noTRs;

                    if (avgL < 25) {

                        System.out.println(r.getPdbCode() + "_" + r.getPdbChain());

                    }
                }

            } catch (Exception ex) {

            }


        }


    }

    public void showFrancois() throws Exception {
        // read data
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");
        for (Row row : rowsFrancois) {
            if (row.getByEyeTot() == 2)
                System.out.println(row.getProtein());

        }

    }


    public void evaluate2() throws Exception {
        List<String> lines = DataIO.readLines("data/evaluation/francois.o");
        Map<String, List<Data>> map = new HashMap<String, List<Data>>();
        for (String line : lines) {
            //System.out.println(line);
            if (line.length() < 10)
                continue;


            String[] data = line.split("\t");
            int start = 0;
            int end = 0;
            if (data[5].startsWith("-")) {
                data[5] = data[5].substring(1, data[5].length());
                start = (-1) * Integer.parseInt(data[5].split("-")[0]);
                end = Integer.parseInt(data[5].split("-")[1]);

            } else {
                start = Integer.parseInt(data[5].split("-")[0]);
                end = Integer.parseInt(data[5].split("-")[1]);
            }


            try {
                //System.out.println(line);
                Data dataobj = new Data();
                dataobj.setPdbBegin(start);
                dataobj.setPdbEnd(end);
                dataobj.setPdbEntry(data[0]);
                dataobj.setScore1(Double.parseDouble(data[7]));
                dataobj.setScore2(Double.parseDouble(data[8]));


                if (map.containsKey(data[0])) {
                    List<Data> l = map.get(dataobj.getPdbEntry());
                    l.add(dataobj);
                    map.put(dataobj.getPdbEntry(), l);
                } else {
                    List<Data> l = new ArrayList<Data>();
                    l.add(dataobj);
                    map.put(dataobj.getPdbEntry(), l);
                }

            } catch (Exception ex) {
                //System.out.println();
            }
        }


        // read from francois
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");
        for (Row row : rowsFrancois) {

            if (row.getByEyeTot() == 1) {

                if (map.containsKey(row.getProtein().trim())) {
                    List<Data> l = map.get(row.getProtein().trim());

                    double maxQA = -1;
                    double maxScore1 = -1;
                    double maxScore2 = -1;
                    for (Data data : l) {
                        maxQA = Math.max(maxQA, data.getQAScore());
                        maxScore1 = Math.max(maxScore1, data.getScore1());
                        maxScore2 = Math.max(maxScore2, data.getScore2());
                    }

                    System.out.println(row.getProtein() + "\t1\t1\t" + NumberFormatUtils.format(maxQA));
//                    boolean check = false;
//                    List<Data> list = map.get(row.getProtein());
//
//                    for (Data data : list) {
//                        int start = data.getPdbBegin();
//                        int end = data.getPdbEnd();
//                        int match = 0;
//                        for (int i = row.getPdbBegin(); i <= row.getPdbEnd(); i++) {
//                            if (start < i && i < end)
//                                match++;
//                        }
//                        double percent = (double) match / (row.getPdbEnd() - row.getPdbBegin() + 1);
//                        if (percent > 0.5) {
//                            check = true;
//                            break;
//                        }
//                    }
//
//                    if (check)
//                        System.out.println(row.getProtein() + "\t1\t1");
//                    else
//                        System.out.println(row.getProtein() + "\t1\t0");
                } else {

                    System.out.println(row.getProtein() + "\t1\t0\t0");

                }

            } else if (row.getByEyeTot() == 2) {

                if (map.containsKey(row.getProtein())) {
                    List<Data> l = map.get(row.getProtein().trim());

                    double maxQA = -1;
                    double maxScore1 = -1;
                    double maxScore2 = -1;
                    for (Data data : l) {
                        maxQA = Math.max(maxQA, data.getQAScore());
                        maxScore1 = Math.max(maxScore1, data.getScore1());
                        maxScore2 = Math.max(maxScore2, data.getScore2());
                    }
                    System.out.println(row.getProtein() + "\t0\t1\t" + NumberFormatUtils.format(maxQA));
                } else {
                    System.out.println(row.getProtein() + "\t0\t0\t0");
                }

            }

        }


    }

    public void evaluate() throws Exception {

        PrintWriter writer = new PrintWriter("data/evaluation.o");
        writer.write("protein\tfinder\tgroup\tnTRs\tavgL\tLR\tRUs\tQA\tR\tclass\tprecise\tconclusion\n");
        Map<String, Row> map = this.readDataFromFrancois();
        //read data from tapo
        List<String> lines = DataIO.readLines("data/francois.o");
        for (String line : lines) {
            //System.out.println(line);
            if (line.length() < 10)
                continue;
            //System.out.println(line);
            String[] data = line.split("\t");
            // System.out.println(data[0]);
            //System.out.println(data[5]);
            // data

            try {
                if (map.containsKey(data[0])) {
                    Row row = map.get(data[0]);
                    int start = 0;
                    int end = 0;
                    if (data[5].startsWith("-")) {
                        data[5] = data[5].substring(1, data[5].length());
                        start = (-1) * Integer.parseInt(data[5].split("-")[0]);
                        end = Integer.parseInt(data[5].split("-")[1]);

                    } else {
                        start = Integer.parseInt(data[5].split("-")[0]);
                        end = Integer.parseInt(data[5].split("-")[1]);
                    }
                    int match = 0;
                    for (int i = row.getPdbBegin(); i <= row.getPdbEnd(); i++) {
                        if (start < i && i < end)
                            match++;
                    }
                    double percent = (double) match / (row.getPdbEnd() - row.getPdbBegin() + 1);

                    if (row.getByEyeTot() == 2)
                        writer.write(line + "\tnoTRs\t0.0\tF\n");
                    else {

                        if (percent > 0.20)
                            writer.write(line + "\tTRs\t" + NumberFormatUtils.format(percent) + "\tT\n");
                        else
                            writer.write(line + "\tTRs\t" + NumberFormatUtils.format(percent) + "\tF\n");
                    }
                }
            } catch (Exception ex) {
                System.out.println(line);
            }

        }


        writer.close();
    }

    public Map<String, Row> readDataFromFrancois() throws Exception {

        // read data

        ProteinCSVReader csvReader = new ProteinCSVReader();

        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");
        Map<String, Row> map = new HashMap<String, Row>();

        for (Row row : rowsFrancois) {

            if (row.getByEyeTot() == 1 || row.getByEyeTot() == 2)
                map.put(row.getProtein(), row);

        }

        return map;
    }


}
