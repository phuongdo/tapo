package org.cnrs.crbm.lib.analysis;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 11/7/2014.
 */
public class NonTRAnalysis {

    public static void main(String[] args) throws Exception {
//        new NonTRAnalysis().analysis();

        //new NonTRAnalysis().extractData();

        //  new NonTRAnalysis().analysisTMScoreForShort2TRs();
//        new NonTRAnalysis().analysisQAForShortTRs();

        ProteinCSVReader csvReader = new ProteinCSVReader();

        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        List<Row> rows = new ArrayList<Row>();
        //convert to rows
        for (RowRepeatDB r : rowsRB) {


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

                    RepeatFinder finder = new RepeatFinder(r.getPdbCode(), r.getPdbChain());
                    // extract repeat units
                    Atom[] atoms = finder.getAtoms();
                    String strSS = finder.getStrSS();
                    List<Repeat> allTRs = new ArrayList<Repeat>();

                    for (int nTRs = 2; nTRs <= 5; nTRs++) {
                        for (int winsize = 10; winsize < 25; winsize = winsize + 1) {
                            for (int i = 0; i < atoms.length - (winsize * nTRs); i = i + 5) {
                                Repeat repeat = new Repeat();
                                // build a segment with at least 2 secondary structure.
                                for (int n = 0; n < nTRs; n++) {
                                    int start = i + n * winsize;
                                    int end = i + (n + 1) * winsize - 1;
                                    String pattern1 = VectorShape.getSSPattern(strSS
                                            .substring(start, end));
                                    if (pattern1.length() >= 2)
                                        repeat.getRepeats().add(new RepeatContent(start, end));
                                }
                                if (repeat.getRepeats().size() >= 2)
                                    allTRs.add(repeat);

                            }
                        }
                    }


                }
            }
        }
    }

    public void analysisTMScoreForShort2TRs() throws Exception {
        // read from francois
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");

        PrintWriter writer = new PrintWriter("data/train/threshold_TRsSmaller25ResAfter+1000.train");
        writer.write("protein\tclass\trmsd\ttmScore\tpattern\tavgL\n");
        MutilAlign align = new MutilAlign();


        for (Row row : rowsFrancois) {
            if (row.getByEyeTot() == 2) {


                int classTRs = row.getByEyeTot();
                if (classTRs == 2)
                    classTRs = 0;
                RepeatFinder finder = new RepeatFinder(row.getPdbCode(), row.getPdbChain());
                // extract repeat units
                Atom[] atoms = finder.getAtoms();
                String strSS = finder.getStrSS();

                double maxTMscore = -1.0;
                double rmsd = 0.0;
                String pattern = "";
                int leng = 0;

                for (int t = 10; t < 25; t = t + 1) {
                    for (int i = 0; i < atoms.length - 2 * t; i = i + 1) {
                        try {
                            Atom[] s1 = Fragement.getFragementsofAtoms(atoms, i, i + t - 1);
                            Atom[] s2 = Fragement.getFragementsofAtoms(atoms, i + t, i + 2 * t - 1);


                            String pattern1 = VectorShape.getSSPattern(strSS
                                    .substring(i, i + t - 1));
                            String pattern2 = VectorShape.getSSPattern(strSS
                                    .substring(i + t, i + 2 * t - 1));

                            int minL = Math.min(s1.length, s2.length);
                            if (pattern1.length() >= 2 && pattern2.equals(pattern1)) {
                                AFPChain afpChain = align.pairAlign(s1, s2);
                                if (afpChain.getOptLength() > 0) {
                                    if (maxTMscore < afpChain.getTMScore()) {
                                        maxTMscore = afpChain.getTMScore();
                                        rmsd = afpChain.getChainRmsd();
                                        pattern = pattern1;
                                        leng = minL;
                                    }

//                                    String segment
//                                            = row.getProtein() + "(" + s1[0].getGroup().getResidueNumber() + "-" + s1[s1.length - 1].getGroup().getResidueNumber() + ":" + s2[0].getGroup().getResidueNumber() + "-" + s2[s2.length - 1].getGroup().getResidueNumber() + ")";
//                                    // System.out.println(minL + "\t" + afpChain.getOptLength() + "\t" + afpChain.getChainRmsd() + "\t" + afpChain.getTMScore());
                                    // process here?
                                    //String content = segment + "\t" + minL + "\t" + NumberFormatUtils.format(afpChain.getChainRmsd()) + "\t" + NumberFormatUtils.format(afpChain.getTMScore()) + "\n";
                                    //writer.write(content);

                                }
                            }


                        } catch (Exception ex) {
                            ex.printStackTrace();

                        }


                    }

                }

                if (maxTMscore > 0) {
                    String content = row.getProtein() + "\t" + classTRs + "\t" + NumberFormatUtils.format(rmsd) + "\t" + NumberFormatUtils.format(maxTMscore) + "\t" + pattern + "\t" + leng +
                            "\n";
                    writer.write(content);
                }

            }


        }


        // for repeats


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
                        RepeatFinder finder = new RepeatFinder(r.getPdbCode(), r.getPdbChain());
                        // extract repeat units
                        Atom[] atoms = finder.getAtoms();
                        String strSS = finder.getStrSS();

                        double maxTMscore = 0.0;
                        double rmsd = 0.0;
                        String pattern = "";

                        int leng = 0;
                        //process
                        //System.out.println(r);
                        for (int t = 10; t < 25; t = t + 1) {
                            for (int i = 0; i < atoms.length - 2 * t; i = i + 1) {
                                try {
                                    Atom[] s1 = Fragement.getFragementsofAtoms(atoms, i, i + t - 1);
                                    Atom[] s2 = Fragement.getFragementsofAtoms(atoms, i + t, i + 2 * t - 1);

                                    String pattern1 = VectorShape.getSSPattern(strSS
                                            .substring(i, i + t - 1));
                                    String pattern2 = VectorShape.getSSPattern(strSS
                                            .substring(i + t, i + 2 * t - 1));

                                    int minL = Math.min(s1.length, s2.length);
                                    if (pattern1.length() >= 2 && pattern2.equals(pattern1)) {
                                        AFPChain afpChain = align.pairAlign(s1, s2);
                                        if (afpChain.getOptLength() > 0) {
                                            if (maxTMscore < afpChain.getTMScore()) {
                                                maxTMscore = afpChain.getTMScore();
                                                pattern = pattern1;
                                                leng = minL;
                                            }

//                                    String segment
//                                            = row.getProtein() + "(" + s1[0].getGroup().getResidueNumber() + "-" + s1[s1.length - 1].getGroup().getResidueNumber() + ":" + s2[0].getGroup().getResidueNumber() + "-" + s2[s2.length - 1].getGroup().getResidueNumber() + ")";
//                                    // System.out.println(minL + "\t" + afpChain.getOptLength() + "\t" + afpChain.getChainRmsd() + "\t" + afpChain.getTMScore());
                                            // process here?
                                            //String content = segment + "\t" + minL + "\t" + NumberFormatUtils.format(afpChain.getChainRmsd()) + "\t" + NumberFormatUtils.format(afpChain.getTMScore()) + "\n";
                                            //writer.write(content);

                                        }
                                    }


                                } catch (Exception ex) {
                                    ex.printStackTrace();

                                }


                            }

                        }
                        if (maxTMscore > 0) {
                            String content = r.getEntry() + "\t1\t" + NumberFormatUtils.format(rmsd) + "\t" + NumberFormatUtils.format(maxTMscore) + "\t" + pattern + "\t" + leng + "\n";
                            writer.write(content);
                        }


                    }
                }

            } catch (Exception ex) {

            }


        }


        writer.close();

    }

    public void analysisQAForShortTRs() throws Exception {
        // read from francois
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");
        for (Row row : rowsFrancois) {
            if (row.getByEyeTot() == 2) {
                RepeatFinder finder = new RepeatFinder(row.getPdbCode(), row.getPdbChain());
                // extract repeat units
                Atom[] atoms = finder.getAtoms();
                String strSS = finder.getStrSS();
                List<Repeat> allTRs = new ArrayList<Repeat>();

                for (int nTRs = 3; nTRs <= 4; nTRs++) {
                    for (int winsize = 10; winsize < 25; winsize = winsize + 1) {
                        for (int i = 0; i < atoms.length - (winsize * nTRs); i = i + 5) {
                            Repeat repeat = new Repeat();
                            // build a segment with at least 2 secondary structure.
                            for (int n = 0; n < nTRs; n++) {
                                int start = i + n * winsize;
                                int end = i + (n + 1) * winsize - 1;
                                String pattern1 = VectorShape.getSSPattern(strSS
                                        .substring(start, end));
                                if (pattern1.length() >= 2)
                                    repeat.getRepeats().add(new RepeatContent(start, end));
                            }
                            if (repeat.getRepeats().size() >= 2)
                                allTRs.add(repeat);

                        }
                    }
                }


                try {

                    double bestQA = -1;
                    double bestQAOld = -1;
                    double bestScore1 = -1;
                    double bestScore2 = -1;
                    double bestScore3 = -1;
                    double bestScore4 = -1;

                    double bestRank = -1;
                    double bestRank2 = -1;

                    for (Repeat repeat : allTRs) {
                        StringBuffer buffer = new StringBuffer();
                        for (RepeatContent content : repeat.getRepeats()) {
                            buffer.append(content.getStart() + "-" + content.getEnd() + ";");
                        }

                        String align = buffer.toString();
                        align = align.substring(0, align.length() - 1);

                        MutilAlign mutilAlign = new MutilAlign(finder.getFeatures(), align);
                        repeat.setScore(mutilAlign.getQAScore());

                        repeat.setScore1(mutilAlign.getScore1());
                        repeat.setScore2(mutilAlign.getScore2());
                        double b = 1;
                        repeat.setRankScore(repeat.getScore() * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));

                        bestQA = Math.max(repeat.getScore(), bestQA);
                        bestScore1 = Math.max(repeat.getScore1(), bestScore1);
                        bestScore2 = Math.max(repeat.getScore2(), bestScore2);
                        bestScore3 = Math.max(mutilAlign.getScore3(), bestScore3);
                        bestScore4 = Math.max(mutilAlign.getScore4(), bestScore4);

                        bestRank = Math.max(repeat.getRankScore(), bestRank);
                        bestQAOld = Math.max(((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3), bestQAOld);
                        bestRank2 = Math.max(((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3) * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b), bestRank2)
                        ;
                        //repeat.setRankScore(repeat.getScore() * Math.pow(1 - 1 / (1 + repeat.getRepeats().size()), a) * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));
                        //System.out.println(bestScore2);
                    }

                    System.out.println(finder.getPdbCode() + "_" + finder.getPdbChain() + "\t" + 0 + "\t" + NumberFormatUtils.format(bestQAOld) + "\t" + NumberFormatUtils.format(bestQA) + "\t" + NumberFormatUtils.format(bestScore1) + "\t" + NumberFormatUtils.format(bestScore2) + "\t" + NumberFormatUtils.format(bestScore3) + "\t" + NumberFormatUtils.format(bestScore4) + "\t" + NumberFormatUtils.format(bestRank) + "\t" + NumberFormatUtils.format(bestRank2));
                } catch (Exception ex) {

                }


            }


        }


        // for repeats


        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        List<Row> rows = new ArrayList<Row>();
        //convert to rows
        for (RowRepeatDB r : rowsRB) {


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

                    RepeatFinder finder = new RepeatFinder(r.getPdbCode(), r.getPdbChain());
                    // extract repeat units
                    Atom[] atoms = finder.getAtoms();
                    String strSS = finder.getStrSS();
                    List<Repeat> allTRs = new ArrayList<Repeat>();

                    for (int nTRs = 2; nTRs <= 5; nTRs++) {
                        for (int winsize = 10; winsize < 25; winsize = winsize + 1) {
                            for (int i = 0; i < atoms.length - (winsize * nTRs); i = i + 5) {
                                Repeat repeat = new Repeat();
                                // build a segment with at least 2 secondary structure.
                                for (int n = 0; n < nTRs; n++) {
                                    int start = i + n * winsize;
                                    int end = i + (n + 1) * winsize - 1;
                                    String pattern1 = VectorShape.getSSPattern(strSS
                                            .substring(start, end));
                                    if (pattern1.length() >= 2)
                                        repeat.getRepeats().add(new RepeatContent(start, end));
                                }
                                if (repeat.getRepeats().size() >= 2)
                                    allTRs.add(repeat);

                            }
                        }
                    }


                    try {

                        double bestQA = -1;
                        double bestQAOld = -1;
                        double bestScore1 = -1;
                        double bestScore2 = -1;
                        double bestScore3 = -1;
                        double bestScore4 = -1;

                        double bestRank = -1;
                        double bestRank2 = -1;

                        for (Repeat repeat : allTRs) {
                            StringBuffer buffer = new StringBuffer();
                            for (RepeatContent content : repeat.getRepeats()) {
                                buffer.append(content.getStart() + "-" + content.getEnd() + ";");
                            }

                            String align = buffer.toString();
                            align = align.substring(0, align.length() - 1);

                            MutilAlign mutilAlign = new MutilAlign(finder.getFeatures(), align);
                            repeat.setScore(mutilAlign.getQAScore());

                            repeat.setScore1(mutilAlign.getScore1());
                            repeat.setScore2(mutilAlign.getScore2());
                            double b = 1;
                            repeat.setRankScore(repeat.getScore() * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));

                            bestQA = Math.max(repeat.getScore(), bestQA);
                            bestScore1 = Math.max(repeat.getScore1(), bestScore1);
                            bestScore2 = Math.max(repeat.getScore2(), bestScore2);
                            bestScore3 = Math.max(mutilAlign.getScore3(), bestScore3);
                            bestScore4 = Math.max(mutilAlign.getScore4(), bestScore4);

                            bestRank = Math.max(repeat.getRankScore(), bestRank);
                            bestQAOld = Math.max(((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3), bestQAOld);
                            bestRank2 = Math.max(((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3) * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b), bestRank2)
                            ;
                            //repeat.setRankScore(repeat.getScore() * Math.pow(1 - 1 / (1 + repeat.getRepeats().size()), a) * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));
                            //System.out.println(bestScore2);
                        }

                        System.out.println(finder.getPdbCode() + "_" + finder.getPdbChain() + "\t" + 1 + "\t" + NumberFormatUtils.format(bestQAOld) + "\t" + NumberFormatUtils.format(bestQA) + "\t" + NumberFormatUtils.format(bestScore1) + "\t" + NumberFormatUtils.format(bestScore2) + "\t" + NumberFormatUtils.format(bestScore3) + "\t" + NumberFormatUtils.format(bestScore4) + "\t" + NumberFormatUtils.format(bestRank) + "\t" + NumberFormatUtils.format(bestRank2));
                    } catch (Exception ex) {

                    }
                }
            }

        }


    }


    public void extractData() throws Exception {
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rows = csvReader.getData("data/RepeatDatalastest.csv");
        List<String> ignoreCase = new ArrayList<String>();
        ignoreCase.add("1yo8");
        ignoreCase.add("2qcs");
        ignoreCase.add("4f5v");
        ignoreCase.add("3m7h");

        for (Row row : rows) {

            if (row.getByEyeTot() != 2)
                continue;

            if (ignoreCase.contains(row.getPdbCode()))
                continue;

            System.out.println(row.getProtein());

        }

    }


    /**
     * We take 2 continuous fragments, sliding and superimposing together from proteins which are
     * non TRs. The size of fragments comes from 10,11, to 90. We also store some importance variable obtained from
     * those superimpose algorithm.
     * With protein with TRs, the same procedure also was applied. Then, we display 2 datasets in the same figure,
     * one with blue plot, others with red plot. This result will confirm our conclusion about what is the best threshold that
     * will be used to identify the TRs to avoid the similar structure happen by accident.
     *
     * @throws Exception
     */
    public void analysis() throws Exception {

        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rows = csvReader.getData("data/RepeatDatalastest.csv");
//        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");
//
//        List<Row> rows = new ArrayList<Row>();
//        //convert to rows
//        for (RowRepeatDB r : rowsRB) {
//
//            if (r.getAnnlevel().equals("Detailed")) {
//                Row row = new Row();
//                row.setProtein(r.getPdbCode() + "_" + r.getPdbChain());
//                rows.add(row);
//            }
//        }

        System.out.println(rows.size());
        //System.exit(1);

        MutilAlign align = new MutilAlign();
        Superimposer superimposer = new Superimposer();


        PrintWriter writer = new PrintWriter("data/threshold_TRs.train");
        writer.write("segment\ta\tg\taRMSD\ttmScore\tnScore\trmsd\tclass\n");
        ProgressBar bar = new ProgressBar();

        List<String> ignoreCase = new ArrayList<String>();
        ignoreCase.add("1yo8");
        ignoreCase.add("2qcs");
        ignoreCase.add("4f5v");
        ignoreCase.add("3m7h");
        bar.update(0, rows.size());
        int process = 1;
        for (Row row : rows) {
            bar.update(process, rows.size());
            process++;

            if (row.getByEyeTot() == 2 || row.getByEyeTot() == 1) {

                int classTR = 1;
                if (row.getByEyeTot() == 2)
                    classTR = 0;

//            if (row.getByEyeTot() != 2)
//                continue;
//
//            if (ignoreCase.contains(row.getPdbCode()))
//                continue;
//            if (!row.getProtein().equals("1isu_A"))
//                continue;
                // System.out.println(row);


                RepeatFinder repeatFinder = new RepeatFinder(row.getPdbCode(), row.getPdbChain());
                // extract repeat units
                Atom[] atoms = repeatFinder.getAtoms();
                String strSS = repeatFinder.getStrSS();


                List<Atom> lAtoms = new ArrayList<Atom>();

                for (Atom atom : atoms)
                    lAtoms.add(atom);
                double minRMSD = Double.MAX_VALUE;
                double minGaps = Double.MAX_VALUE;
                double maxScore = Double.MIN_VALUE;
                double maxTM = Double.MIN_VALUE;
                int maxOptLength = Integer.MIN_VALUE;

                double minRMSD_Real = Double.MAX_VALUE;

                for (int t = 10; t < 100; t = t + 1) {

                    //System.out.println(t);
                    /**
                     * for each of window size, we take only the minimum RMSD
                     */

                    for (int i = 0; i < atoms.length - 2 * t; i = i + 1) {

                        try {
                            Atom[] s1 = Fragement.getFragementsofAtoms(atoms, i, i + t - 1);
                            Atom[] s2 = Fragement.getFragementsofAtoms(atoms, i + t, i + 2 * t - 1);

                            String pattern1 = VectorShape.getSSPattern(strSS
                                    .substring(i, i + t - 1));
                            String pattern2 = VectorShape.getSSPattern(strSS
                                    .substring(i + t, i + 2 * t - 1));

                            if (pattern1.length() >= 2 && pattern1.equals(pattern2)) {

                                int minL = Math.min(s1.length, s2.length);
                                int maxL = Math.max(s1.length, s2.length);
                                //System.out.println(minL + " and " + maxL);
                                AFPChain afpChain = align.pairAlign(s1, s2);

                                if (afpChain.getOptLength() > 0) {
                                    String segment
                                            = row.getProtein() + "(" + s1[0].getGroup().getResidueNumber() + "-" + s1[s1.length - 1].getGroup().getResidueNumber() + ":" + s2[0].getGroup().getResidueNumber() + "-" + s2[s2.length - 1].getGroup().getResidueNumber() + ")";
                                    // System.out.println(minL + "\t" + afpChain.getOptLength() + "\t" + afpChain.getChainRmsd() + "\t" + afpChain.getTMScore());
                                    // process here?

                                    minRMSD = Math.min(afpChain.getChainRmsd(), minRMSD);
                                    maxOptLength = Math.max(afpChain.getOptLength(), maxOptLength);
                                    minGaps = Math.min(afpChain.getGapLen(), minGaps);
                                    maxTM = Math.max(afpChain.getTMScore(), maxTM);
                                    maxScore = Math.max(afpChain.getAlignScore(), maxScore);
                                    //System.out.println(superimposer.superimposeSimple(s1, s2));
                                    minRMSD_Real = Math.min(superimposer.superimposeSimple(s1, s2), minRMSD_Real);
                                    //String content = segment + "\t" + minL + "\t" + afpChain.getOptLength() + "\t" + afpChain.getGapLen() + "\t" + NumberFormatUtils.format(afpChain.getChainRmsd()) + "\t" + NumberFormatUtils.format(afpChain.getTMScore()) + "\n";

                                }
                            }
                        } catch (Exception ex) {
                            ex.printStackTrace();

                        }


                    }


//
//                List<List<Atom>> subSets = Lists.partition(lAtoms,
//                        t);

                    //System.out.println(subSets.size());
                    // int start = this.getPosition(atoms, row.getPdbBegin());
                    // int end = this.getPosition(atoms, row.getPdbEnd());


                    // int leng = (int) Math.floor((double) (end - start + 1) / row.getN());

//
//                List<Atom[]> units = new ArrayList<Atom[]>();
//                for (List<Atom> latom1 : subSets) {
//
//                    Atom[] atomSet1 = latom1.toArray(new Atom[latom1.size()]);
//                    units.add(atomSet1);
//                }

//                for (int i = start; i < end && i + leng < end; i = i + leng) {
//
//                    //System.out.println(i + "-" + (i + leng - 1));
//                    Atom[] unit = Fragement.getFragementsofAtoms(atoms, i, i + leng - 1);
//                    units.add(unit);
//                }


//                for (int i = 0; i < units.size(); i++) {
//                    for (int j = i + 1; j < units.size(); j++) {
//
//                        Atom[] s1 = units.get(i);
//                        Atom[] s2 = units.get(j);
//
//                            int minL = Math.min(s1.length, s2.length);
//                            int maxL = Math.max(s1.length, s2.length);
//                            AFPChain afpChain = pairAlign.pairAlign(s1, s2);
//
//                            if (afpChain.getOptLength() > 0) {
//
//
//                                String segment
//                                        = row.getProtein() + "(" + s1[0].getGroup().getResidueNumber() + "-" + s1[s1.length - 1].getGroup().getResidueNumber() + ":" + s2[0].getGroup().getResidueNumber() + "-" + s2[s2.length - 1].getGroup().getResidueNumber() + ")";
//                                // System.out.println(minL + "\t" + afpChain.getOptLength() + "\t" + afpChain.getChainRmsd() + "\t" + afpChain.getTMScore());
//                                String content = segment + "\t" + minL + "\t" + afpChain.getOptLength() + "\t" + afpChain.getGapLen() + "\t" + NumberFormatUtils.format(afpChain.getChainRmsd()) + "\t" + NumberFormatUtils.format(afpChain.getTMScore()) + "\n";
//
//                                if (content.length() > 10)
//                                    writer.write(content);
//                            }
//                        } catch (Exception ex) {
//
//                        }
//
//                    }
//
//                }


                }

                // write here
                String segment
                        = row.getProtein();
                String content = segment + "\t" + maxOptLength + "\t" + minGaps + "\t" + NumberFormatUtils.format(minRMSD) + "\t" + NumberFormatUtils.format(maxTM) + "\t" + NumberFormatUtils.format(maxScore) + "\t" + NumberFormatUtils.format(minRMSD_Real) + "\t" + classTR + "\n";
                if (content.length() > 10 && content.length() < 50) {
                    //System.out.println(content);
                    writer.write(content);
                }


            }
        }

        writer.close();


    }


    int getPosition(Atom[] atoms, int posNsq) {

        for (int i = 0; i < atoms.length; i++) {

            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;

        }

        return 0;
    }

}
