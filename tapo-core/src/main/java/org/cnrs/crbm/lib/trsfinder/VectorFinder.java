package org.cnrs.crbm.lib.trsfinder;

import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.conf.ThresholdConfig;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;
import org.cnrs.crbm.lib.repeats.module.VectorModule;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

public class VectorFinder extends Finder {

    public VectorFinder(Features features) {
        super(features);
        this.name = "VECTORS";

    }

    final static int MIN_ELEMENTS = 2;
    final static int MAX_ELEMENTS = 8;

    double maxVScore = 0.0;
    static Logger logger = LoggerFactory.getLogger(VectorFinder.class);

    @Override
    public void findRepeat(Features features) {
        try {

//            double maxTMScore = 0.0;
            //VectorModule vectorModule = new VectorModule();
            List<ProVector> ssVectors = ProVector.toSecondaryVector(features
                    .getVectors());

            TMEvaluation tmEvaluation = new TMEvaluation();
            for (int winsize = MIN_ELEMENTS; winsize <= MAX_ELEMENTS; winsize++) {

                List<Repeat> lstRepeat = this.getRepeats(ssVectors, features.getAtoms(),
                        winsize);

                repeats.addAll(lstRepeat);
            }
//            this.getCombineScore().setTmScore(maxTMScore);
            this.getCombineScore().setvScore(maxVScore);

        } catch (Exception e) {
            e.printStackTrace();
            logger.error(e.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }

        // return repeats;

    }


    public List<Repeat> getRepeats(List<ProVector> vectors, Atom[] atoms,
                                   int winsize) throws Exception {
        List<Repeat> newRepeats = new ArrayList<Repeat>();
        Superimposer superimpose = new Superimposer();

        List<List<ProVector>> list = new ArrayList<List<ProVector>>();
        for (int i = 0; i < vectors.size() - winsize + 1; i++) {
            List<ProVector> sub_list = new ArrayList<ProVector>();
            for (int j = 0; j < winsize; j++) {
                ProVector v = vectors.get(i + j);
                sub_list.add(v);
            }
            list.add(sub_list);
        }
        int lastTRs = -1;
        for (int j = 0; j < list.size(); j++) {

            if (lastTRs >= 0 && lastTRs + winsize > j) {
                j = lastTRs + winsize - 1;
                continue;
            }
            List<ProVector> sub_list_ref = list.get(j);
            String pattern_ref = this.getPatterOfVectors(sub_list_ref);


            int start1 = sub_list_ref.get(0).getPosStart();
            int end1 = sub_list_ref.get(sub_list_ref.size() - 1).getPosEnd();

            Repeat repeat = new Repeat();
            repeat.getRepeats().add(new RepeatContent(start1, end1));
            /**
             *scan to the right
             */
            for (int i = j + winsize; i < list.size(); ) {
                List<ProVector> sub_list = list.get(i);
                // check pattern
                String pattern_sub = this.getPatterOfVectors(sub_list);
                if (pattern_sub.equals(pattern_ref)) {
                    double vscore = superimpose.compareListVector(sub_list_ref,
                            sub_list);
                    //check only 2 adjacent set of vectors
                    if (i == (j + winsize)) {
                        maxVScore = Math.max(maxVScore, vscore);
                    }
                    if (vscore > ThresholdConfig.VECTOR_THRES) {
                        repeat.getRepeats().add(
                                new RepeatContent(sub_list.get(0).getPosStart(),
                                        sub_list.get(sub_list.size() - 1)
                                                .getPosEnd()));
                        i = i + winsize;
                        lastTRs = i;
                    } else {
                        //i = i + 1;
                        break;
                    }
                } else {
                    //i = i + 1;
                    break;
                }
            }// end scan to the right

            /**
             *scan to the left
             */
            for (int i = j - winsize; i > 0; ) {
                List<ProVector> sub_list = list.get(i);
                // check pattern
                String pattern_sub = this.getPatterOfVectors(sub_list);
                if (pattern_sub.equals(pattern_ref)) {
                    double vscore = superimpose.compareListVector(sub_list_ref,
                            sub_list);
                    //check only 2 adjacent set of vectors
                    if (i == (j - winsize)) {
                        maxVScore = Math.max(maxVScore, vscore);
                    }
                    if (vscore > ThresholdConfig.VECTOR_THRES) {
                        repeat.getRepeats().add(
                                new RepeatContent(sub_list.get(0).getPosStart(),
                                        sub_list.get(sub_list.size() - 1)
                                                .getPosEnd()));
                        i = i - winsize;
                    } else {
                        //i = i - 1;
                        break;

                    }
                } else {
                    //i = i - 1;
                    break;
                }
            }// end scan to the left


            if (repeat.getRepeats().size() >= 2) {
                repeat.sortByPosition();
                newRepeats.add(repeat);

            }

        }


        return newRepeats;

    }

    private String getPatterOfVectors(List<ProVector> sub_list) {
        String pattern_sub = "";
        for (ProVector v : sub_list) {
            pattern_sub += v.getType();
        }
        return pattern_sub;
    }
}
