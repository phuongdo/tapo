package org.cnrs.crbm.lib.statistic;

/**
 * Created by pdoviet on 1/28/2015.
 */
public class Data {

    private String pdbEntry;
    private int pdbBegin;
    private int pdbEnd;
    private double score1 = 0.0;
    private double score2 = 0.0;
    private int nTRs = 0;

    public double getScore1() {
        return score1;
    }

    public void setScore1(double score1) {
        this.score1 = score1;
    }

    public int getnTRs() {
        return nTRs;
    }

    public void setnTRs(int nTRs) {
        this.nTRs = nTRs;
    }

    public double getScore2() {
        return score2;
    }

    public void setScore2(double score2) {
        this.score2 = score2;
    }

    public double getQAScore() {

        return (score1 + score2) / 2;
    }

    public String getPdbEntry() {
        return pdbEntry;
    }

    public void setPdbEntry(String pdbEntry) {
        this.pdbEntry = pdbEntry;
    }

    public int getPdbBegin() {
        return pdbBegin;
    }

    public void setPdbBegin(int pdbBegin) {
        this.pdbBegin = pdbBegin;
    }

    public int getPdbEnd() {
        return pdbEnd;
    }

    public void setPdbEnd(int pdbEnd) {
        this.pdbEnd = pdbEnd;
    }
}
