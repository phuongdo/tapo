package org.cnrs.crbm.lib.io;

import org.cnrs.crbm.lib.repeats.CombineScore;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 7/17/2015.
 */
public class TaPoFastaFormat {

    private String header = "";
    private boolean is3DRepeat = false;
    private String pdbID = "";
    private CombineScore combineScore;

    public double getSvmScore() {
        return svmScore;
    }

    private double svmScore;
    private List<Repeat> repeats = new ArrayList<Repeat>();


    public TaPoFastaFormat(String aheader, String adescription) {
        // convert to its own format
        this.header = aheader;

        // parse header
        this.pdbID = aheader.split("\\|")[0];
        this.combineScore = new CombineScore(aheader.split("\\|")[1]);
        this.svmScore = Double.parseDouble(aheader.split("\\|")[2].split("=")[1]);
        // String adescription
        if (!adescription.contains("No-TRs")) {
            this.is3DRepeat = true;
            String[] lines = adescription.split("\\|");
            for (String line : lines) {
                if (line.trim().length() > 0) {
                    Repeat repeat = this.convertToRepeat(line);
                    this.repeats.add(repeat);
                }
            }

        }
    }

    public Repeat convertToRepeat(String line) {
        Repeat repeat = new Repeat();
        String[] data = line.split("\t");
        try {

            String[] units = data[6].split(";");
            try {
                repeat.setScore(Double.parseDouble(data[7]));
                repeat.setRankScore(Double.parseDouble(data[8]));
            } catch (Exception ex) {
            }
            repeat.setFinderName(data[1]);
            repeat.setCluster(data[2]);


            for (String unit : units) {
                int start = 0;
                int end = 0;
                if (unit.startsWith("-")) {
                    unit = unit.substring(1, unit.length());
                    start = (-1) * Integer.parseInt(unit.split("-")[0]);
                    end = Integer.parseInt(unit.split("-")[1]);
                    repeat.getRepeats().add(new RepeatContent(start, end));

                } else {
                    start = Integer.parseInt(unit.split("-")[0]);
                    end = Integer.parseInt(unit.split("-")[1]);
                    repeat.getRepeats().add(new RepeatContent(start, end));
                }
            }

        } catch (Exception ex) {

        }

        return repeat;

    }

    public String getPdbID() {
        return pdbID;
    }

    public void setPdbID(String pdbID) {
        this.pdbID = pdbID;
    }

    public CombineScore getCombineScore() {
        return combineScore;
    }

    public void setCombineScore(CombineScore combineScore) {
        this.combineScore = combineScore;
    }

    public String getHeader() {
        return header;
    }

    public void setHeader(String header) {
        this.header = header;
    }

    public List<Repeat> getRepeats() {
        return repeats;
    }

    public void setRepeats(List<Repeat> repeats) {
        this.repeats = repeats;
    }

    public boolean is3DRepeat() {
        return is3DRepeat;
    }

    public void setIs3DRepeat(boolean is3DRepeat) {
        this.is3DRepeat = is3DRepeat;
    }
}
