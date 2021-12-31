package org.cnrs.crbm.lib.repeats.clusters;

import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Created by pdoviet on 10/10/2014.
 */
public class ClusterRepeat {


    public List<Repeat> getRepeats() {
        return repeats;
    }

    public void setRepeats(List<Repeat> repeats) {
        this.repeats = repeats;
    }

    private List<Repeat> repeats = new ArrayList<Repeat>();

    public Region getConcensusRegion() {
        return concensusRegion;
    }

    public void setConcensusRegion(Region concensusRegion) {
        this.concensusRegion = concensusRegion;
    }

    private Region concensusRegion;

    public Repeat getLongestRepeat() {

        Repeat repeatMax = null;

        int max = -1;
        for (Repeat repeat : repeats) {

            int size = repeat.getEnd() - repeat.getStart();
            if (max < size) {
                max = size;
                repeatMax = repeat;
            }
        }
        return repeatMax;

    }

    public void assignRScore() {
        int minStart = 100;
        int maxEnd = 0;
        int maxNoTRs = 0;
        double minRU = 1000;
        for (Repeat repeat : this.repeats) {
            int start = repeat.getStart();
            int end = repeat.getEnd();
            minStart = Math.min(start, minStart);
            maxEnd = Math.max(end, maxEnd);
            maxNoTRs = Math.max(maxNoTRs, repeat.getRepeats().size());
            //minRU = Math.min(minRU, repeat.getAvgLength());

        }
        //int size = Math.abs(maxEnd - minStart);
        int expectedSize = this.concensusRegion.size();
        int expectedNoTRs = maxNoTRs;// (int) (size / minRU);
        for (Repeat repeat : this.repeats) {
            int start = repeat.getStart();
            int end = repeat.getEnd();
            repeat.setRankScore(repeat.getScore() * ((double) (end - start + 1) / expectedSize) * ((double) repeat.getRepeats().size() / expectedNoTRs));
            //System.out.println(repeat.getRankScore());
        }

    }


    public int size() {
        return repeats.size();
    }

    public void sortScoreDESC() {
        Collections.sort(repeats, new Comparator<Repeat>() {
            public int compare(Repeat s1, Repeat s2) {
                if (s1.getRankScore() > s2.getRankScore())
                    return -1;
                else if (s1.getRankScore() < s2.getRankScore())
                    return 1;
                else return 0;
            }
        });

    }

    @Override
    public String toString() {
        StringBuffer buffer = new StringBuffer();
        for (Repeat repeat : repeats) {
            buffer.append(repeat.getStart() + "-" + repeat.getEnd() + " " + repeat.getScore() + "\n");
        }
        return buffer.toString();
    }
}
