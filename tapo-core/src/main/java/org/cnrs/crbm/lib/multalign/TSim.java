package org.cnrs.crbm.lib.multalign;

import org.cnrs.crbm.lib.math.MapUtil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 10/9/2014.
 */
public class TSim {


    public static void main(String[] args) {

//        List<String> mas = new ArrayList<String>();
//        mas.add("ETCLL-DILDTAG-Q-------E-EYSAMRDQYMRT-----");
//        mas.add("------GFLCVFAIN-------NTKSFEDIHQYREQIKRVK");
//        mas.add("---DVPMVLVGNKC-DLAGRTV----ESRQAQDLARSY---");
//        mas.add("------IPYIETSA-KT-----R-QGVEDAFYTLVREIR--");

        List<String> msa = new ArrayList<String>();
        msa.add("0010000111111111000----000000");
        msa.add("00000111110000111000000001000");
        msa.add("00000001110---111000000100000");
        TSim tsim = new TSim();

        System.err.println(tsim.getHETATMScore(msa));


    }


    public double getTsimScore(List<String> mas) {

        int size = mas.get(0).length();
        if (size == 0)
            return 0.0;
        List<Map<Character, Double>> profile = new ArrayList<Map<Character, Double>>();
        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < size; i++) {
            Map<Character, Double> freqs = new HashMap<Character, Double>();

            for (int j = 0; j < mas.size(); j++) {
                Character strClass = mas.get(j).charAt(i);
                if (freqs.containsKey(strClass)) {
                    freqs.put(strClass, freqs.get(strClass) + 1.00);
                } else {
                    freqs.put(strClass, 1.00);
                }
            }
            for (Map.Entry<Character, Double> entry : freqs.entrySet()) {
                entry.setValue(entry.getValue() / mas.size());
            }
            profile.add(freqs);
            freqs = MapUtil.sortByValue(freqs);
            buffer.append(freqs.keySet().toArray()[0]);
            //String consensus = (String) freqs.keySet().toArray()[0];
        }

        // compare profile with msa
        String consensus = buffer.toString();

        double tsim = 0.0;
        double D = 0.0;
        double N = size * mas.size();
        for (String seq : mas) {

            //System.out.println("> string x ");
            // compare and get score
            double score = 0.0;
            for (int i = 0; i < seq.length(); i++) {

                char c = seq.charAt(i);
                Map<Character, Double> f = profile.get(i);
                if (c != '-')
                    score += f.get(c);
            }
            score = score / seq.length();
            //System.out.println(score);
            D += score;
            //D+= Calcs.hammingDistance(consensus,seq);
        }

        if (mas.size() > 0)
            tsim = D / mas.size();

        return tsim;

    }


    /**
     * calculation of the consensus region of HETATM ION score
     * 0010000111111111000----000000
     * 00000111110000111000000001000
     * 00000001110---111000000100000
     *
     * @param mas
     * @return
     */
    public double getHETATMScore(List<String> mas) {
        int size = mas.get(0).length();
        if (size == 0)
            return 0.0;
        double tsim = 0.0;
        // find the largest but > 0

        int[][] matrix = new int[mas.size()][size];

        int nrows = mas.size();
        int ncols = size;
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                if (mas.get(i).charAt(j) != '-') {
                    matrix[i][j] = Integer.parseInt("" + mas.get(i).charAt(j));
                } else
                    matrix[i][j] = -1;
            }
        }


        int maxH = 0;
        double totalPercentageCore = 0;
        for (int i = 0; i < ncols; i++) {
            int colsScore = 0;
            for (int j = 0; j < nrows; j++) {
                if (matrix[j][i] != -1) {
                    colsScore += matrix[j][i];
                }
            }
//            if ((double) colsScore / nrows > 0.4) {
//                totalPercentageCore++;
//            }

            double consensus = (double) colsScore / nrows;
            totalPercentageCore+=consensus;
            if (consensus> 0.0) {
                maxH++;
            }

        }


//        for (int i = 0; i < nrows ; i++) {
//            int total=0;
//            for (int j = 0; j < ncols ; j++) {
//                if(matrix[i][j]==1){
//                    total++;
//                }
//            }
//
//            maxH = Math.max(maxH,total);
//        }
        tsim = (double) totalPercentageCore / maxH;
        return tsim;

    }


}
