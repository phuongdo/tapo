package org.cnrs.crbm.lib.repeats;

import org.cnrs.crbm.lib.utils.NumberFormatUtils;

/**
 * Created by pdoviet on 5/30/2015.
 */
public class CombineScore {


    public static void main(String[] args) {
        CombineScore combineScore = new CombineScore("TM-Score=0.419;Psim-Score=1.000;CE-Score=8.341;V-Score=0.241;L-Score=0.000;CA-Score=0.000;S-Score=0.083;RA-Score=0.400");
        System.out.println(combineScore);
    }

    public CombineScore() {

    }

    public static String toHeader() {
        return "Vscore,TMscore,CEscore,CAscore,Lscore,Sscore,Psimscore,RAscore";
    }

    public CombineScore(String strScores) {

//        System.out.println(strScores);
        String[] scores = strScores.split(";");
        this.tmScore = Double.parseDouble(scores[0].split("=")[1]);
        this.psimScore = Double.parseDouble(scores[1].split("=")[1]);
        this.ceScore = Double.parseDouble(scores[2].split("=")[1]);
        this.vScore = Double.parseDouble(scores[3].split("=")[1]);
        this.lScore = Double.parseDouble(scores[4].split("=")[1]);
        this.caScore = Double.parseDouble(scores[5].split("=")[1]);
        try {
            this.sigScore = Double.parseDouble(scores[6].split("=")[1]);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
        this.rapScore = Double.parseDouble(scores[7].split("=")[1]);
    }

    public double getTmScore() {
        return tmScore;
    }

    public void setTmScore(double tmScore) {
        this.tmScore = tmScore;
    }

    public double getPsimScore() {
        return psimScore;
    }

    public void setPsimScore(double psimScore) {
        this.psimScore = psimScore;
    }

    public double getCeScore() {
        return ceScore;
    }

    public void setCeScore(double ceScore) {
        this.ceScore = ceScore;
    }

    public double getvScore() {
        return vScore;
    }

    public void setvScore(double vScore) {
        this.vScore = vScore;
    }

    public double getRapScore() {
        return rapScore;
    }

    public void setRapScore(double rapScore) {
        this.rapScore = rapScore;
    }

    private double ceScore = 0.0; // CE-Symm Score
    private double vScore = 0.0; // V-Score
    private double tmScore = 0.0;// max TM-Score
    private double psimScore = 0.0;// max Conformation Alphabet Score
    private double caScore = 0.0; // {0,1}; matched score form predefined patterns
    private double lScore = 0.0; // H-Score from hetam module
    private double sigScore = 0.0; // did not support now


    private double rapScore = 0.0;//


    public double getCaScore() {
        return caScore;
    }

    public void setCaScore(double caScore) {
        this.caScore = caScore;
    }

    public double getlScore() {
        return lScore;
    }

    public void setlScore(double lScore) {
        this.lScore = lScore;
    }

    public double getSigScore() {
        return sigScore;
    }

    public void setSigScore(double sigScore) {
        this.sigScore = sigScore;
    }


    @Override
    public String toString() {
        return "TM-Score=" + NumberFormatUtils.format(this.tmScore)
                + ";Psim-Score=" + NumberFormatUtils.format(this.psimScore)
                + ";CE-Score=" + NumberFormatUtils.format(this.ceScore)
                + ";V-Score=" + NumberFormatUtils.format(this.vScore)
                + ";L-Score=" + NumberFormatUtils.format(this.lScore)
                + ";CA-Score=" + NumberFormatUtils.format(this.caScore)
                + ";S-Score=" + NumberFormatUtils.format(this.sigScore)
                + ";RA-Score=" + NumberFormatUtils.format(this.rapScore);
    }

    public String toCSVFormat() {
        return NumberFormatUtils.format(this.vScore)
                + "," + NumberFormatUtils.format(this.tmScore)
                + "," + NumberFormatUtils.format(this.ceScore)
                + "," + NumberFormatUtils.format(this.caScore)
                + "," + NumberFormatUtils.format(this.lScore)
                + "," + NumberFormatUtils.format(this.sigScore)
                + "," + NumberFormatUtils.format(this.psimScore)
                + "," + NumberFormatUtils.format(this.rapScore);
    }

    public String toLibSVMFormat() {
        return "1 "
                + "1:" + NumberFormatUtils.format(this.vScore)
                + " 2:" + NumberFormatUtils.format(this.tmScore)
                + " 3:" + NumberFormatUtils.format(this.ceScore)
                + " 4:" + NumberFormatUtils.format(this.caScore)
                + " 5:" + NumberFormatUtils.format(this.lScore)
                + " 6:" + NumberFormatUtils.format(this.sigScore)
                + " 7:" + NumberFormatUtils.format(this.psimScore)
                + " 8:" + NumberFormatUtils.format(this.rapScore);
    }

    public String toLibSVMFormat(String label) {
        return label + " "
                + "1:" + NumberFormatUtils.format(this.vScore)
                + " 2:" + NumberFormatUtils.format(this.tmScore)
                + " 3:" + NumberFormatUtils.format(this.ceScore)
                + " 4:" + NumberFormatUtils.format(this.caScore)
                + " 5:" + NumberFormatUtils.format(this.lScore)
                + " 6:" + NumberFormatUtils.format(this.sigScore)
                + " 7:" + NumberFormatUtils.format(this.psimScore)
                + " 8:" + NumberFormatUtils.format(this.rapScore);
    }


}

