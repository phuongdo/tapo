package org.cnrs.crbm.lib.analysis;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 11/5/2014.
 */
public class SubMatrixAnalysis {


    public static void main(String[] args) {

        SubMatrixAnalysis analysis = new SubMatrixAnalysis();
        analysis.readSubMatrix();


    }


    void readSubMatrix() {

        // read the files
        final File folder = new File(Dir.MUSTANG_TMP_DIR + "/msa_output/");
        Map<String, Integer> pairs = listFilesForFolder(folder);
        String letters = "EKLIBVPSDCGHAFXY".toLowerCase();
        //String letters = "ailstv".toLowerCase();

        // cal Fkm

        int Fkm = 0;

//        for (Map.Entry<String, Integer> entry : pairs.entrySet()) {
//            Fkm += entry.getValue();
//        }
//
//        System.out.println(Fkm);
        int a = 0;
        double[][] Q = new double[letters.length()][letters.length()];
        double[] P = new double[letters.length()];
        double[][] E = new double[letters.length()][letters.length()];
        double[][] S = new double[letters.length()][letters.length()];
        double[][] F = new double[letters.length()][letters.length()];
        for (int i = 0; i < letters.length(); i++) {
            for (int j = i; j < letters.length(); j++) {
                char c1 = letters.charAt(i);
                char c2 = letters.charAt(j);
                String key = "";
                if (c1 < c2)
                    key = c1 + "" + c2;
                else
                    key = c2 + "" + c1;
                int freq = 0;
                if (pairs.containsKey(key))
                    freq = pairs.get(key);
                //System.out.println(key + ":" + freq);

                //double Qij = (double) freq / Fkm;
                //System.out.println(Qij);
                //Q[i][j] = Qij;
                F[i][j] = freq;
                F[j][i] = freq;


            }

        }


        for (int i = 0; i < F.length; i++) {
            for (int j = 0; j < i; j++) {
                Fkm += F[i][j];

            }
        }
        Fkm = Fkm * 2;
        //System.out.println(Fkm);

        for (int i = 0; i < Q.length; i++) {
            for (int j = i; j < Q.length; j++) {
                Q[i][j] = (double) F[i][j] / Fkm;
                Q[j][i] = Q[i][j];

            }
        }


        // System.out.println(a);

        for (int i = 0; i < P.length; i++) {
            double tmp = 0.0;
            for (int j = 0; j < P.length; j++) {
                if (i != j) {
                    tmp += Q[i][j];

                }

            }


            P[i] = Q[i][i] + tmp / 2;

        }


        for (int i = 0; i < letters.length(); i++) {
            for (int j = i; j < letters.length(); j++) {

                if (i == j) {
                    E[i][j] = P[i] * P[j];
                } else {
                    E[i][j] = 2 * P[i] * P[j];
                    E[j][i] = E[i][j];
                }
            }

        }

        for (int i = 0; i < letters.length(); i++) {
            for (int j = i; j < letters.length(); j++) {

                //System.out.println(Q[i][j] + ":" + E[i][j]);
                if (Q[i][j] == 0 || E[i][j] == 0)
                    S[i][j] = 0;
                else
                    S[i][j] = Math.round(2 * log2(Q[i][j] / E[i][j]));
                S[j][i] = S[i][j];
            }

        }


        for (int i = 0; i < S.length; i++) {
            for (int j = 0; j < S.length; j++) {

//                if (S[i][j] != 0)
//                    S[i][j] = S[i][j] + 20;
//                else
//                    S[i][j] = 0;
                System.out.print((int) S[i][j] + "\t");
            }
            System.out.println();
        }

    }

    static double log2(double x) {
        return (Math.log10(x) / Math.log10(2));
    }

    public Map<String, Integer> listFilesForFolder(final File folder) {

        Map<String, Integer> pairs = new HashMap<String, Integer>();

        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isDirectory()) {
                listFilesForFolder(fileEntry);
            } else {
                //System.out.println(fileEntry.getAbsolutePath());
                List<String> msa = DataIO.readLines(fileEntry.getAbsolutePath());

                int nRows = msa.size();
                int nCols = msa.get(0).length();

                for (int i = 0; i < nCols; i++) {
                    for (int j = 0; j < nRows; j++) {
                        for (int k = j + 1; k < nRows; k++) {
                            char c1 = msa.get(j).charAt(i);
                            char c2 = msa.get(k).charAt(i);
                            String key = "";
                            if (c1 < c2)
                                key = c1 + "" + c2;
                            else
                                key = c2 + "" + c1;

                            if (!pairs.containsKey(key))
                                pairs.put(key, 1);
                            else {
                                pairs.put(key, pairs.get(key) + 1);
                            }

                        }

                    }

                }

            }
        }

        return pairs;
    }

}
