package org.cnrs.crbm.lib.repeats;

import org.cnrs.crbm.lib.trsfinder.Region;

/**
 * Created by pdoviet on 8/28/2015.
 */
public class RepeatRegion extends Region {

    private int predRU;

    public int getPredRU() {
        return predRU;
    }

    public void setPredRU(int predRU) {
        this.predRU = predRU;
    }

    public int getPredLen() {
        return predLen;
    }

    public void setPredLen(int predLen) {
        this.predLen = predLen;
    }

    private int predLen;
    public RepeatRegion(int start, int end) {
        super(start, end);
    }

    public int getSeqStart() {
        return seqStart;
    }

    public void setSeqStart(int seqStart) {
        this.seqStart = seqStart;
    }

    public int getSeqEnd() {
        return seqEnd;
    }

    public void setSeqEnd(int seqEnd) {
        this.seqEnd = seqEnd;
    }

    private int seqStart;
    private int seqEnd;



}
