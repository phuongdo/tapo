package org.cnrs.crbm.lib.repeats.clusters;

import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 9/10/2015.
 */
public class ConcensusRegion {
    final static double THRESHOLD_CUTOFF = 0.3;

    /**
     * @param lstRep list of repeats
     * @param L      length of the protein structure
     * @return
     */

    public static List<Region> findConcensusRegion(List<Repeat> lstRep, int L) {
        double[] consensus = new double[L];
        double[] consensusSize = new double[L];
        for (int i = 0; i < consensusSize.length; i++) {
            int size = 0;
            for (Repeat repeat : lstRep) {
                if (repeat.getStart() < i && i < repeat.getEnd()) {
                    size++;
                }
            }
            consensusSize[i] = size;
        }
        for (Repeat repeat : lstRep) {
            for (int i = repeat.getStart(); i < repeat.getEnd(); i++) {
                consensus[i] += repeat.getScore();
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
        return findRegion(consensus);

    }


    public static List<Region> findRegion(double[] concencus) {
        List<Region> lstRegions = new ArrayList<Region>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < concencus.length; i++) {
            double ch = concencus[i];
            if (ch >= THRESHOLD_CUTOFF) {
                start = i;
                while (i < concencus.length
                        && concencus[i] >= THRESHOLD_CUTOFF) {
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

}
