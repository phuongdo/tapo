package org.cnrs.crbm.lib.cm;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.*;
import java.util.List;


public class ContactMap {

    public String getSecondaryStructure() {
        return secondaryStructure;
    }

    public void setSecondaryStructure(String secondaryStructure) {
        this.secondaryStructure = secondaryStructure;
    }

    private String pdbCode = "";
    private double cutoff = 7.00;
    private String CT = "Cb";
    private int noRes = 0;
    private Atom[] setAtom = null;
    private double[] histogram = null;
    private String secondaryStructure = "";
    private List<Pair> cmaps = new ArrayList<Pair>();

    public Map<Integer, List<Integer>> getCmapsList() {
        return cmapsList;
    }

    private Map<Integer, List<Integer>> cmapsList = new HashMap<Integer, List<Integer>>();

    //private int[][] cmapMatrix = null;

    public void setCmaps(List<Pair> cmaps) {
        this.cmaps = cmaps;
    }

    public ContactMap(Atom[] setAtom, double cutoff, String secondaryStructure) {

        this.setAtom = setAtom;
        this.cutoff = cutoff;
        noRes = this.setAtom.length;
        histogram = new double[noRes];

        //cmapMatrix = new int[noRes][noRes];
        this.secondaryStructure = secondaryStructure;
        this.calContactMap();
    }

    public ContactMap(Atom[] setAtom, double cutoff) {
        this.setAtom = setAtom;
        this.cutoff = cutoff;
        noRes = this.setAtom.length;
        histogram = new double[noRes];
        this.calContactMap();
    }


    private void calContactMap() {
        // cal contact Map here

        for (int i = 0; i < noRes - 1; i++) {
            for (int j = i; j < noRes; j++) {

                Atom atom1 = setAtom[i];
                Atom atom2 = setAtom[j];

                // System.out.println("i:" + i + "j:" + j + "dis"
                // + Calc.getDistance(atom1, atom2));
                if (Calc.getDistance(atom1, atom2) <= this.cutoff) {
                    Pair pair = new Pair();
                    pair.setFirst(i);
                    pair.setSecond(j);
                    cmaps.add(pair);

                    if (cmapsList.containsKey(i)) {
                        List<Integer> l = cmapsList.get(i);
                        l.add(j);
                        cmapsList.put(i, l);
                    } else {
                        List<Integer> l = new ArrayList<Integer>();
                        l.add(j);
                        cmapsList.put(i, l);
                    }

                    if (cmapsList.containsKey(j)) {
                        List<Integer> l = cmapsList.get(j);
                        l.add(i);
                        cmapsList.put(j, l);
                    } else {
                        List<Integer> l = new ArrayList<Integer>();
                        l.add(i);
                        cmapsList.put(j, l);
                    }

                    // cal distance
                }

            }
        }
    }

    /**
     * We are interested in contact map overlap to find the repetitive because
     * it tell us how the pattern of one protein is folding - Contact Map
     * CM(i,j) < Delta ( with Delta < 8A) - IF secondary structure is +HELIX :
     * only consider i and j with (j >=i+5) +BETA STRAND: i and j with (j >=i
     * +4) IF TURN do nothing Update: (28-03) - in case of turn, we do not
     * consider to the contact map - HHHHH <=5 and BB were consider it as turn.
     */
    private void calContactMapWithSecondaryStructure() {
        // cal contact Map here

        for (int i = 0; i < noRes - 1; i++) {

            char ss_i = secondaryStructure.charAt(i);
            if (ss_i == '-')
                continue;
            int distance = 0;
            if (ss_i == 'H')
                distance = 5;
            else if (ss_i == 'B')
                distance = 4;
            else {
                distance = 1;
            }
            for (int j = i + distance; j < noRes; j++) {
                char ss_j = secondaryStructure.charAt(j);
                if (ss_j == '-')
                    continue;
                Atom atom1 = setAtom[i];
                Atom atom2 = setAtom[j];

                try {
                    // System.out.println("i:" + i + "j:" + j + "dis"
                    // + Calc.getDistance(atom1, atom2));
                    if (Calc.getDistance(atom1, atom2) <= this.cutoff) {
                        Pair pair = new Pair();
                        pair.setFirst(i);
                        pair.setSecond(j);
                        cmaps.add(pair);

                        if (cmapsList.containsKey(i)) {
                            List<Integer> l = cmapsList.get(i);
                            l.add(j);
                            cmapsList.put(i, l);
                        } else {
                            List<Integer> l = new ArrayList<Integer>();
                            l.add(j);
                            cmapsList.put(i, l);
                        }

                        if (cmapsList.containsKey(j)) {
                            List<Integer> l = cmapsList.get(j);
                            l.add(i);
                            cmapsList.put(j, l);
                        } else {
                            List<Integer> l = new ArrayList<Integer>();
                            l.add(i);
                            cmapsList.put(j, l);
                        }

                        // cal distance
                    }

                } catch (Exception e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
        }
    }

    public String toFile() {
        StringBuilder builder = new StringBuilder();
        builder.append("#PDB: " + this.pdbCode + "\n");
        builder.append("#CT: " + this.CT + "\n");
        builder.append("#CUTOFF: " + this.cutoff + "\n");
        for (Pair pair : cmaps) {
            builder.append(pair.getFirst() + " " + pair.getSecond() + "\n");
        }
        return builder.toString();

    }

    public double[] getHistogram() {
        int length = setAtom.length;
        for (int i = 0; i < length; i++) {
            if (cmapsList.containsKey(i)) {
                List<Integer> list = cmapsList.get(i);
                for (Integer j : list) {
                    if (j > i) {
                        int delta = j - i + 1;
                        histogram[i] = Math.max(histogram[i], delta);
                    }
                }

            }
        }
        return histogram;

    }

    public double[] getHistogramWithRange(int max) {
        //int length = setAtom.length;

        for (Pair pair : cmaps) {
            int i = pair.getFirst();
            int j = pair.getSecond();
            int delta = j - i + 1;
            if (max > delta) {
                histogram[i] = Math.max(histogram[i], delta);
                histogram[j] = Math.max(histogram[j], delta);
            }
        }

        return histogram;
    }

    public List<Pair> getCmaps() {
        return cmaps;
    }

    public int getNoRes() {
        return this.noRes;
    }

    public String getPdbCode() {
        return pdbCode;
    }

    public void setPdbCode(String pdbCode) {
        this.pdbCode = pdbCode;
    }

}
