package org.cnrs.crbm.lib.seqalign;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 5/23/2015.
 */
public class RepeatSeq {

    public List<RepeatUnitSeq> getUnits() {
        return units;
    }
    public void setUnits(List<RepeatUnitSeq> units) {
        this.units = units;
    }
    List<RepeatUnitSeq>  units = new ArrayList<RepeatUnitSeq>();

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    double score =0.0; // if exists


}
