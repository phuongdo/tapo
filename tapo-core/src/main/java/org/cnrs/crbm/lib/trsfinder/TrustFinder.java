package org.cnrs.crbm.lib.trsfinder;

import org.cnrs.crbm.lib.multalign.PSim;
import org.cnrs.crbm.lib.multalign.TSim;
import org.cnrs.crbm.lib.repeats.shortTRs.ShortTRsFilter;
import org.cnrs.crbm.lib.seqalign.RepeatSeq;
import org.cnrs.crbm.lib.seqalign.RepeatUnitSeq;
import org.cnrs.crbm.lib.seqalign.Trust;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

public class TrustFinder extends Finder {

    public TrustFinder(Features features) {
        super(features);
        this.name = "TRUST-23A";

    }

    static Logger logger = LoggerFactory.getLogger(TrustFinder.class);

    @Override
    public void findRepeat(Features features) {

        try {
            double maxPSim = 0.0;
            String pattern = this.features.getStrSeqSADB();
            int gapo = -8;
            int gape = -2;
            Trust trust = new Trust("conf/BLOSUM62",
                    "ABCDEFGHIKLMNPQRSTVWXYZO", this.features.getPdbCode(),
                    pattern, gapo, gape);
            TSim tSim = new TSim();
            List<RepeatSeq> repeatsSq = trust.getRepeats();
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
                double score = pSim.getSimilarity();

                // check repeat
                a_repeat.sortByPosition();
                if (ShortTRsFilter.filter(a_repeat, features.getStrSeqAp())) {
                    maxPSim = Math.max(score, maxPSim);
                    if (score >= 0.7) {
                        this.repeats.add(a_repeat);
                    }
                }

            }

            // store maxPsim
            this.getCombineScore().setPsimScore(maxPSim);


        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());

        }

    }
}
