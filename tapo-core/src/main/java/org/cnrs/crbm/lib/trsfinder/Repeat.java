package org.cnrs.crbm.lib.trsfinder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * repeat object
 *
 * @author phuongdo
 */
public class Repeat {


    private double score1;

    public double getScore2() {
        return score2;
    }

    public void setScore2(double score2) {
        this.score2 = score2;
    }

    public double getScore1() {
        return score1;
    }

    public void setScore1(double score1) {
        this.score1 = score1;
    }

    private double score2;

    public double getScore3() {
        return score3;
    }

    public void setScore3(double score3) {
        this.score3 = score3;
    }

    private double score3;

    public double getScore4() {
        return score4;
    }

    public void setScore4(double score4) {
        this.score4 = score4;
    }

    private double score4;

    public String getCluster() {
        return cluster;
    }

    public void setCluster(String cluster) {
        this.cluster = cluster;
    }

    private String cluster;

    List<RepeatContent> repeats = new ArrayList<RepeatContent>();

    public List<RepeatContent> getRepeats() {
        return repeats;
    }

    public void setRepeats(List<RepeatContent> repeats) {
        this.repeats = repeats;
    }

    public double getAvgLength() {
        double avg = 0.0;
        for (RepeatContent r : repeats) {
            avg += r.size();
        }
        return avg / repeats.size();
    }

    public int getStart() {
        return this.repeats.get(0).getStart();
    }

    public int getEnd() {
        return this.repeats.get(this.repeats.size() - 1).getEnd();

    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    private double score = 0.0;
    private double rankScore = 0.0;

    public String getFinderName() {
        return finderName;
    }

    public void setFinderName(String finderName) {
        this.finderName = finderName;
    }

    private String finderName;

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        for (RepeatContent rc : repeats)
            builder.append(rc.toString() + "\n");
        return builder.toString();
    }

    public static void sortByLengthDES(List<Repeat> rs) {

        Collections.sort(rs, new Comparator<Repeat>() {
            public int compare(Repeat s1, Repeat s2) {
                if (s1.getAvgLength() > s2.getAvgLength())
                    return -1;
                else
                    return 1;
            }
        });

    }

    public void sortByPosition() {
        Collections.sort(repeats, new Comparator<RepeatContent>() {
            public int compare(RepeatContent s1, RepeatContent s2) {
                if (s1.getStart() > s2.getStart())
                    return 1;
                else
                    return -1;
            }
        });

    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Repeat) {
            Repeat other = (Repeat) obj;
            if (other.getRepeats().size() == this.getRepeats().size() && other.getStart() == this.getStart() && other.getEnd() == this.getEnd())
                return true;
            else return false;

        } else {
            return false;
        }
    }


    public String getRepeatRegionString() {

        StringBuffer buffer = new StringBuffer();
        double avgL = 0.0;
        for (RepeatContent rc : repeats) {
            buffer.append(rc.getStart() + "-" + rc.getEnd() + ";");
        }
        String frags = buffer.toString();
        frags = frags.substring(0, frags.length() - 1);

        return frags;
    }

    public double getRankScore() {
        return rankScore;
    }

    public void setRankScore(double rankScore) {
        this.rankScore = rankScore;
    }
}
