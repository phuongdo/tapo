package org.cnrs.crbm.lib.trsfinder;

import org.cnrs.crbm.lib.multalign.PSim;
import org.cnrs.crbm.lib.repeats.shortTRs.ShortTRsFilter;
import org.cnrs.crbm.lib.seqalign.RepeatSeq;
import org.cnrs.crbm.lib.seqalign.RepeatUnitSeq;
import org.cnrs.crbm.lib.seqalign.TREksWapper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;


public class TReksSADBFinder extends Finder {
    static Logger logger = LoggerFactory.getLogger(TReksSADBFinder.class);

    public TReksSADBFinder(Features features) {
        super(features);
        this.name = "TREKS-23A";
    }

    @Override
    public void findRepeat(Features features) {
        String pattern = features.getStrSeqSADB();
        try {
            double maxPSim = 0.0;
            TREksWapper TReks = new TREksWapper("", pattern);
            List<RepeatSeq> repeatsSq = TReks.getRepeatSeqs();
            for (RepeatSeq repeat : repeatsSq) {
                // each of repeats
                //abc
                Repeat a_repeat = new Repeat();
                List<String> msa = new ArrayList<String>();
                for (RepeatUnitSeq ru : repeat.getUnits()) {
                    msa.add(ru.getStrAlign());
                    a_repeat.getRepeats().add(new RepeatContent(ru.getStart(), ru.getEnd()));
                }
                PSim pSim = new PSim(msa);
                a_repeat.sortByPosition();
                if (ShortTRsFilter.filter(a_repeat, features.getStrSeqAp()) && a_repeat.getAvgLength() >= 10) {
                    // check repeat
                    maxPSim = Math.max(maxPSim, pSim.getSimilarity());
                    if (pSim.getSimilarity() >= 0.7) {
                        this.repeats.add(a_repeat);
                    }
                }

            }

            this.combineScore.setPsimScore(maxPSim);

        } catch (Exception ex) {
            // ex.printStackTrace();
            logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }

        //this.repeats = repeats_copy;

    }
}
