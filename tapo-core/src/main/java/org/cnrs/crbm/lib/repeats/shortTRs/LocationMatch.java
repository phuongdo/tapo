package org.cnrs.crbm.lib.repeats.shortTRs;

/**
 * Created by pdoviet on 1/26/2015.
 */
public class LocationMatch {

    private int start;
    private int end;
    private String meg = "";


    public int getStart() {
        return start;
    }

    public LocationMatch(int start, int end, String meg) {
        this.start = start;
        this.end = end;
        this.meg = meg;
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

    public String getMeg() {
        return meg;
    }

    public void setMeg(String meg) {
        this.meg = meg;
    }
}
