package org.cnrs.crbm.lib.trsfinder;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.raphael.Raphael;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.module.FinderOutput;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;


public class RaphaelFinder extends Finder {

    static Logger logger = LoggerFactory.getLogger(RaphaelFinder.class);

    public RaphaelFinder(Features features) {
        super(features);
        this.name = "RapBase";
//		this.name = "RMSD";

    }

    public void setRegionScan(Region regionScan) {
        this.regionScan = regionScan;
    }

    Region regionScan;

    @Override
    public void findRepeat(Features features) {
        try {
            TMEvaluation tmEvaluation = new TMEvaluation();
            Atom[] atomScan = Fragement.getFragementsofAtoms(atoms, regionScan.getStart(), regionScan.getEnd());
            Raphael raphael = new Raphael(atomScan);
            if (raphael.getRepeatLength() > 0 & raphael.getTotalScore() > 0.05) {
                int repeatLeng = (int) raphael.getRepeatLength();
                FinderOutput finderOut = tmEvaluation.findRepeatsBasedOnTMScore(repeatLeng, regionScan.getStart(), regionScan.getEnd(), features);
                List<Repeat> lrepeats = finderOut.getRepeats();
                if (lrepeats.size() > 0)
                    this.repeats.addAll(lrepeats);
                this.combineScore.setRapScore(raphael.getTotalScore());
                this.combineScore.setTmScore(finderOut.getCombineScore().getTmScore());
            }


        } catch (Exception ex) {
//            ex.printStackTrace();
            logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }

    }

}
