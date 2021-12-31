package org.cnrs.crbm.lib.seqalign;

/**
 * Created by pdoviet on 5/23/2015.
 */
public class RepeatUnitSeq {

    String strAlign = "";

    public RepeatUnitSeq(String strAlign, int start, int end) {
        this.strAlign = strAlign;
        this.start = start;
        this.end = end;
    }

    int start = 0;
    int end = 0;

    public String getStrAlign() {
        return strAlign;
    }

    public void setStrAlign(String strAlign) {
        this.strAlign = strAlign;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int size() {
        return end - start + 1;
    }
}
