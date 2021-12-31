package org.cnrs.crbm.trclassification;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.cnrs.crbm.lib.cm.ContactMap;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 9/30/2015.
 * In this class, we try to investigate
 * tandem repeat protein formed "opened" or "closed" structure.
 * Closed structure is defined as a fixed number of repeats due
 * to their circular or ‘‘closed’’ structures as opposed to "opened" structure in
 * which there are now limitation of growth along axis.
 */
public class OpenStructure {

    public static void main(String[] args) throws Exception {

        new OpenStructure().benchmarkRepeatLength();
    }

    double getContactScore(Atom[] atoms, Repeat repeat) {

        double score = 0;
        RepeatContent rcStart = repeat.getRepeats().get(0);
        RepeatContent rcEnd = repeat.getRepeats().get(repeat.getRepeats().size() - 1);
        int noContact = 0;
        for (int i = rcStart.getStart(); i <= rcStart.getEnd(); i++) {
            Atom atomi = atoms[i];
            for (int j = rcEnd.getStart(); j <= rcEnd.getEnd(); j++) {
                Atom atomj = atoms[j];
                double deltaD = Calc.getDistance(atomi, atomj);
                if (deltaD < 7.0) {
                    noContact++;
                    break;
                }
            }
        }
        return (double) noContact / repeat.getAvgLength();
    }


    double globalScore(Atom[] atoms, Map<Integer, List<Integer>> cmapsList) {
        double score = 0.0;
        //ContactMap contactMap = new ContactMap(atoms, 7.0);
        int L = atoms.length;
        int thresholdL = L / 3;
        // Map<Integer, List<Integer>> cmapsList = contactMap.getCmapsList();
        int noResCM = 0;
        for (int i = 0; i < L; i++) {
            int maxDelta = 0;
            if (cmapsList.containsKey(i)) {
                List<Integer> list = cmapsList.get(i);
                for (Integer j : list) {
                    int delta = Math.abs(j - i) + 1;
                    maxDelta = Math.max(delta, maxDelta);
                }
            }
            if (maxDelta >= thresholdL)
                noResCM++;
        }
        score = (double) noResCM / (2 * thresholdL);

        return score;

    }


    public void benchmarkRepeatLength() throws Exception {

        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/benchmarkset/tapov1.1.1.output");
        ClusterLocation clusterLocation = new ClusterLocation();
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        System.out.println("pdb\tlabel\tdiScore\tcmScore\tglobularScore\tcoverRegion");
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
                if (strClass.contains("UA"))
                    continue;
                if (region.startsWith("-")) {
                    region = region.substring(1, region.length());
                    startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                    endRegion = Integer.parseInt(region.split("-")[1]);

                } else {
                    startRegion = Integer.parseInt(region.split("-")[0]);
                    endRegion = Integer.parseInt(region.split("-")[1]);
                }

                Region refRegion = new Region(startRegion, endRegion);
                double avgRefL = Double.parseDouble(data[3]);
                StringBuffer buffer = new StringBuffer();
                String strStart = pdb + "\t" + refRegion + "\t" + avgRefL + "\t";
                String falsePositive = "0.00\t0-0\t0.00";
                if (fastaTRs.containsKey(pdb)) {
                    TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);
                    if (taPoFastaFormat.is3DRepeat()) {
                        List<Repeat> lstRepeat = new ArrayList<Repeat>();
                        for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                            if (repeat.getCluster().contains("selected")) {
                                lstRepeat.add(repeat);
                            }
                        }
                        // get
                        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                        Atom[] atoms = repeatFinder.getAtoms();
                        //repeatFinder.getCmapsList();

                        String display = "";
                        boolean found = false;
                        for (Repeat repeat : lstRepeat) {

                            int start = PdbTools.getResSeq(atoms, repeat.getStart());
                            int end = PdbTools.getResSeq(atoms, repeat.getEnd());
                            Region predRegion = new Region(start, end);
                            double cover = clusterLocation.cover(predRegion, refRegion);
                            if (cover >= 0.5) {
//                                double scoreContact = this.getContactScore(Fragement.getFragementsofAtoms(atoms, repeat.getStart(), repeat.getEnd()));
                                //double score = this.getContactScore(atoms);
                                double globalScore = this.globalScore(atoms, repeatFinder.getCmapsList());
                                double contactScore = this.getContactScore(atoms, repeat);
                                double score = 0;
                                List<RepeatContent> repeatContents = repeat.getRepeats();

                                // distance from start and end region
                                double deltaD = Calc.getDistance(atoms[repeat.getStart()], atoms[repeat.getRepeats().get(repeat.getRepeats().size() - 1).getStart()]);
                                double delta = 0;

                                for (int i = 0; i < repeat.getRepeats().size() - 1; i++) {
                                    RepeatContent rc1 = repeat.getRepeats().get(i);
                                    RepeatContent rc2 = repeat.getRepeats().get(i + 1);
                                    delta += Calc.getDistance(atoms[rc1.getStart()], atoms[rc2.getStart()]);
                                }

                                if (delta > 0)
                                    score = deltaD / delta;

                                double coverRegion = repeat.getEnd() - repeat.getStart() + 1;
                                coverRegion = coverRegion / atoms.length;
                                String outClass = "open";
                                if (strClass.contains("IV."))
                                    outClass = "close";
                                System.out.println(pdb + "\t" + outClass + "\t" + NumberFormatUtils.format(score) + "\t" + NumberFormatUtils.format(contactScore) + "\t" + NumberFormatUtils.format(globalScore) + "\t" + NumberFormatUtils.format(coverRegion));

                            }

                        }

                    } else {
//                        System.out.println(strStart + "0\t" + falsePositive);
                    }


                } else {
//                    System.out.println(strStart + falsePositive);
                }


//                System.out.println(pdb + "\t" + refRegion + "\t" + avgRefL + "\t" + buffer.toString());

            } catch (Exception ex) {
//                ex.printStackTrace();
            }
        }

    }
}
