package org.cnrs.crbm.lib.trsfinder;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.cnrs.crbm.lib.repeats.SeedVector;
import org.cnrs.crbm.lib.repeats.module.FinderOutput;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

/**
 * Created by pdoviet on 3/14/2015.
 */
public class VectorSignalFinder extends Finder {

    static Logger logger = LoggerFactory.getLogger(VectorSignalFinder.class);

    public VectorSignalFinder(Features features) {
        super(features);
        this.name = "VECTORS";

    }


    @Override
    public void findRepeat(Features features) {
        // List<Repeat> repeats = new ArrayList<Repeat>();

        try {

            TMEvaluation tmEvaluation = new TMEvaluation();
            List<ProVector> ssVectors = ProVector.toSecondaryVector(features
                    .getVectors());
            List<List<Integer>> seeds = SeedVector.scanSeed(ssVectors);
            for (List<Integer> seedPos : seeds) {
                DescriptiveStatistics stats = new DescriptiveStatistics();
                double diff = 0.0;
                for (int pos = 0; pos < seedPos.size() - 1; pos++) {
                    stats.addValue(seedPos.get(pos + 1) - seedPos.get(pos));
                }
                double P_avg = stats.getMean();
                double P_sd = stats.getStandardDeviation();

                Repeat candidate = new Repeat();
                for (int i = 0; i < seedPos.size() - 1; i++) {
                    candidate.getRepeats().add(new RepeatContent(seedPos.get(i), seedPos.get(i + 1) - 1));
                }

                int start_ending_Pos = seedPos.get(seedPos.size() - 1);
                //add end fragment
                int end_ending_Pos = start_ending_Pos + (int) P_avg - 1;
                if (end_ending_Pos >= atoms.length) {
                    end_ending_Pos = atoms.length - 1;
                }
                candidate.getRepeats().add(new RepeatContent(start_ending_Pos, end_ending_Pos));
                FinderOutput finderOutput = tmEvaluation.findRepeatsOfOneSegments(candidate, features);
                List<Repeat> listR = finderOutput.getRepeats();
                this.combineScore.setTmScore(finderOutput.getCombineScore().getTmScore());
                if (listR.size() > 0)
                    this.repeats.addAll(listR);

            }


        } catch (Exception e) {
            e.printStackTrace();
            logger.error(e.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }

        // return repeats;

    }


    private String getPatterOfVectors(List<ProVector> sub_list) {
        String pattern_sub = "";
        for (ProVector v : sub_list) {
            pattern_sub += v.getType();
        }
        return pattern_sub;
    }
}
