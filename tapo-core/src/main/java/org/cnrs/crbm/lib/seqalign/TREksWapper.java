package org.cnrs.crbm.lib.seqalign;

import java.util.ArrayList;
import java.util.List;

import nptr.*;

/**
 * <pre>
 * T-Reks is a novel programme that's very sensitive with short repeats.
 *  It was developed by Andrey Kajava'teams 2010
 *  Website: http://bioinfo.montp.cnrs.fr/?Rlib=t-reks/
 * </pre>
 *
 * @author phuongdo
 */
public class TREksWapper {

    List<Repeat> repeats = new ArrayList<Repeat>();

    @Deprecated
    public List<Repeat> getRepeats() {
        return repeats;
    }

    @Deprecated
    public void setRepeats(List<Repeat> repeats) {
        this.repeats = repeats;
    }

    public List<RepeatSeq> getRepeatSeqs() {
        return repeatSeqs;
    }

    List<RepeatSeq> repeatSeqs = new ArrayList<RepeatSeq>();


    public TREksWapper(String dec, String seq) {
        Parameters.setParamDefault("20");
        SeqRepeat seqRep = new SeqRepeat(dec, seq);
        TRExplorer Trex = new TRExplorer(seqRep);
        OverlapManager aCopies = new OverlapManager();
        Trex.run();
        aCopies = Trex.getACopies();
        if (aCopies.size() >= 1) {
            for (int i = 0; i < aCopies.size(); i++) {
                if (((Repeat) aCopies.get(i)).isReal()) {
                    /**
                     * @Deprecated
                     * WILL BE REMOVE IN THE NEXT VERSION OF TAPO
                     */
                    Repeat repeat = (Repeat) aCopies.get(i);
                    repeats.add(repeat);
                    /**
                     *
                     */
                    RepeatSeq repeatSeq = new RepeatSeq();
                    for (int a = 0; a < repeat.getCopies().size(); ++a) {
                        Copy copy = ((Copy) repeat.getCopies().get(a));
                        RepeatUnitSeq ru = new RepeatUnitSeq(copy.getSequence(), copy.getBeginPosition(), copy.getEndPosition());
                        repeatSeq.getUnits().add(ru);
                        repeatSeqs.add(repeatSeq);
                    }
                    //repeatSeq.getUnits().add(re);

                }
            }
            seqRep.clear();
            aCopies.clear();
        } else {
            // System.out.println("repeat not found in sequence "
            // + seqRep.getDesc());
        }

    }
}
