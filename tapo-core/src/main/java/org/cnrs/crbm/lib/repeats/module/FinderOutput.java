package org.cnrs.crbm.lib.repeats.module;

import org.cnrs.crbm.lib.repeats.CombineScore;
import org.cnrs.crbm.lib.trsfinder.Repeat;

import java.util.ArrayList;
import java.util.List;

/**
 * This class contains TR candidates and score
 * Created by pdoviet on 6/12/2015.
 */
public class FinderOutput {

    CombineScore combineScore = new CombineScore();
    List<Repeat> repeats = new ArrayList<Repeat>();

    public CombineScore getCombineScore() {
        return combineScore;
    }

    public void setCombineScore(CombineScore combineScore) {
        this.combineScore = combineScore;
    }

    public List<Repeat> getRepeats() {
        return repeats;
    }

    public void setRepeats(List<Repeat> repeats) {
        this.repeats = repeats;
    }
}
