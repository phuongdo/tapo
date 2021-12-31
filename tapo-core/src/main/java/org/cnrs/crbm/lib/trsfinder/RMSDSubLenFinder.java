package org.cnrs.crbm.lib.trsfinder;

import org.cnrs.crbm.lib.repeats.module.FinderOutput;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

/**
 * @author pdoviet
 *         from TAPO 1.1.2 we will use greedy method
 *         to scan whole structure which is smaller than 300 residues.
 */
public class RMSDSubLenFinder extends Finder {

    final static boolean debug = false;

    static Logger logger = LoggerFactory.getLogger(RMSDSubLenFinder.class);

    public void setListWinsize(List<Integer> listWinsize) {
        this.listWinsize = listWinsize;
    }

    List<Integer> listWinsize = new ArrayList<Integer>();


    public RMSDSubLenFinder(Features features) {
        super(features);
        this.name = "RMSD";

    }

    @Override
    public void findRepeat(Features features) {
        try {
            TMEvaluation tmEvaluation = new TMEvaluation();
            for (int predLen : listWinsize) {
                FinderOutput finderOut = tmEvaluation.findRepeatsBasedOnTMScore(predLen, 0, atoms.length - 1, features);
                List<Repeat> lrepeats = finderOut.getRepeats();
                if (lrepeats.size() > 0)
                    this.repeats.addAll(lrepeats);
                this.combineScore.setTmScore(finderOut.getCombineScore().getTmScore());
            }

        } catch (Exception ex) {
//            ex.printStackTrace();
            logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }

    }

    public static List<Integer> generateLens(int min, int max) {
        List<Integer> lstLens = new ArrayList<Integer>();
        for (int i = min; i < max; i++) {
            lstLens.add(i);
        }
        return lstLens;

    }

}
