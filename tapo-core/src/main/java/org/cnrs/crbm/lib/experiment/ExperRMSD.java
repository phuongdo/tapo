package org.cnrs.crbm.lib.experiment;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.repeats.*;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.repeats.clusters.ClusterRepeat;
import org.cnrs.crbm.lib.trsfinder.AtomFinder;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 1/8/2015.
 * We want to test how many times will be spent. by using only this
 * RMSD method. See Algorithm 1 in our report.
 */
public class ExperRMSD {


    StringBuilder contentBuilder = new StringBuilder();

    Features features = null;

    public ExperRMSD(String pdbCode, String pdbChain) {
        RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
        features = finder.getFeatures();
    }

    public void findRepeat() {
        List<Repeat> allTRs = new ArrayList<Repeat>();

        for (int winsize = 10; winsize < 90; winsize = winsize + 5) {
            //System.out.println(winsize);
            allTRs.addAll(AtomFinder.getRepeat(features, winsize));

        }

        Atom[] atoms = features.getAtoms();
        String pdbCode=  features.getPdbCode();
        String pdbChain = features.getPdbChain();

        // ranking
        try {
            //cluster

            List<ClusterRepeat> clusters = new ClusterLocation().cluster(allTRs);

            int indexG = 1;
            for (ClusterRepeat c : clusters) {
                // System.out.println("-----");
                c.sortScoreDESC();

                for (Repeat top : c.getRepeats()) {
                    int start = top.getStart();
                    int end = top.getEnd();
                    // pdbId_chain Finder avgLeng nUnits RL
                    int seqStart = atoms[start].getGroup()
                            .getResidueNumber().getSeqNum();
                    int seqEnd = atoms[end].getGroup()
                            .getResidueNumber().getSeqNum();

                    StringBuffer unitsBuffer = new StringBuffer();
                    for (RepeatContent unit : top.getRepeats()) {

                        int seqS = atoms[unit.getStart()].getGroup()
                                .getResidueNumber().getSeqNum();
                        int seqE = atoms[unit.getEnd()].getGroup()
                                .getResidueNumber().getSeqNum();

                        unitsBuffer.append(seqS + "-" + seqE + ";");
                    }
                    String unitsStr = unitsBuffer.toString();
                    unitsStr = unitsStr.substring(0, unitsStr.length() - 1);
                    contentBuilder.append(pdbCode + "_" + pdbChain + "\t" + top.getFinderName() + "\t" + "clus" + indexG
                            + "\t"
                            + top.getRepeats().size() + "\t"
                            + NumberFormatUtils.format(top.getAvgLength()) + "\t" + seqStart
                            + "-" + seqEnd + "\t" + unitsStr + "\t" + NumberFormatUtils.format(top.getScore()) + "\n")
                    ;
                }

                indexG++;

            }

        } catch (Exception ex) {

        }


    }

    public String getOutput(){
        return contentBuilder.toString();
    }

}
