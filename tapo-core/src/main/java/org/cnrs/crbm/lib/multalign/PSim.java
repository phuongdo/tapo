package org.cnrs.crbm.lib.multalign;

/**
 * Created by pdoviet on 5/26/2015.
 */

import java.util.ArrayList;
import java.util.List;

import nptr.utils.Hash;

public class PSim {
    //private AlignCopies alignment;
    private List<String> msa = new ArrayList<String>();
    private String consensus;
    private double similarity;


    public PSim(List<String> msa) {
        this.msa = msa;
        this.compute();
    }


    private void compute() throws StringIndexOutOfBoundsException {
        String pattern = "";
        this.consensus = "";
        int countEmptyCol = 0;
        int nCols = this.msa.get(0).length();
        if (this.msa.size() > 0) {
            double nbMatchs = 0.0D;

            try {
                for (int ex = 0; ex < nCols; ++ex) {
                    Hash column = new Hash();

                    int nbResCons;
                    for (nbResCons = 0; nbResCons < this.msa.size(); ++nbResCons) {
                        pattern = this.msa.get(nbResCons);
                        int nb = 1;

                        try {
                            Character e = new Character(pattern.charAt(ex));
                            if (column.containsKey(e)) {
                                nb = Integer.parseInt(column.get(e).toString()) + 1;
                            }

                            column.put(e, Integer.valueOf(nb));
                        } catch (Exception var10) {
//                            System.out.println("could not pairAlign:");
//                            System.out.println(this.alignment);
                        }
                    }

                    column = column.sort(false);
                    nbResCons = Integer.parseInt(column.getFirstValue().toString());
                    if (nbResCons > 1) {
                        if (!column.getFirst().equals(new Character('-'))) {
                            this.consensus = this.consensus + column.getFirst();
                            nbMatchs += (double) nbResCons;
                        } else if (nbResCons != this.msa.size()) {
                            this.consensus = this.consensus + column.getFirst();
                            nbMatchs += (double) nbResCons;
                        } else {
                            ++countEmptyCol;
                            this.consensus = this.consensus + "-";
                        }
                    } else {
                        this.consensus = this.consensus + "X";
                        ++nbMatchs;
                    }
                }
            } catch (StringIndexOutOfBoundsException var11) {
                var11.printStackTrace();
                //System.exit(1);
            }

            this.similarity = nbMatchs / (double) (this.msa.size() * (nCols - countEmptyCol));
        }

    }


    public String getConsensus() {
        return this.consensus;
    }

    public void setConsensus(String consensus) {
        this.consensus = consensus;
    }

    public double getSimilarity() {
        return this.similarity;
    }

    public void setSimilarity(double similarity) {
        this.similarity = similarity;
    }
}