package org.cnrs.crbm.lib.trsfinder;

import java.util.ArrayList;
import java.util.List;

import org.cnrs.crbm.lib.multalign.PSim;
import org.cnrs.crbm.lib.repeats.shortTRs.ShortTRsFilter;
import org.cnrs.crbm.lib.sadb.Conformation;
import org.cnrs.crbm.lib.sadb.Seq3d;
import org.cnrs.crbm.lib.seqalign.RepeatSeq;
import org.cnrs.crbm.lib.seqalign.RepeatUnitSeq;
import org.cnrs.crbm.lib.seqalign.Trust;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TrustCombinedFinder extends Finder {

    public TrustCombinedFinder(Features features) {
        super(features);
        this.name = "TRUST-16A";

    }

    static Logger logger = LoggerFactory.getLogger(TrustCombinedFinder.class);

    @Override
    public void findRepeat(Features features) {

        Seq3d seq3D = new Seq3d(this.features.getStrSeqAp());
        String pattern = seq3D.getSeq();
        List<Conformation> confors = seq3D.getConfors();
        int gapo = -6;
        int gape = -2;
        // if (this.features.getPdbCode().equals("1BIK"))
        // System.out.println();
        Trust trust = new Trust("conf/BLOSUM_COMBINED", "EKLIBVPSDCGHAFXYOWMN",
                this.features.getPdbCode(), pattern, gapo, gape);
        // System.out.println(pattern);
        try {

            double maxPSim = 0.0;
            //TSim tSim = new TSim();
            List<RepeatSeq> repeatsSq = trust.getRepeats();
            for (RepeatSeq repeat : repeatsSq) {
                // each of repeats
                //abc
                Repeat a_repeat = new Repeat();
                List<String> msa = new ArrayList<String>();
                for (RepeatUnitSeq ru : repeat.getUnits()) {
                    msa.add(ru.getStrAlign());
                    int start = ru.getStart();
                    int end = ru.getEnd();
                    int leng = 0;
                    // convert current to repeats
                    // get real position
                    // count in real leng of resi
                    for (int i = start; i <= end; i++) {
                        leng += confors.get(i).getEnd()
                                - confors.get(i).getStart() + 1;
                    }
                    start = confors.get(start).getStart();
                    a_repeat.getRepeats().add(
                            new RepeatContent(start, start + leng - 1));


                }

                PSim pSim = new PSim(msa);
                double score = pSim.getSimilarity();
                maxPSim = Math.max(score, maxPSim);
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

            logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }

    }
}
