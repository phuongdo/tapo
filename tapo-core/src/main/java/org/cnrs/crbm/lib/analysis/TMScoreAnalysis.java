package org.cnrs.crbm.lib.analysis;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Features;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 2/25/2015.
 */
public class TMScoreAnalysis {


    public static void main(String[] args) {
        String pdbCode = "2o6x";
        String pdbChain = "A";
        if (args.length > 0) {
            pdbCode = args[0];
            pdbChain = args[1];
        }
        try {
//            TMScoreAnalysis evaluation = new TMScoreAnalysis(pdbCode, pdbChain);
//            evaluation.findRepeat();
//            System.out.println(evaluation.getOutput());
            //System.exit(0);

        } catch (Exception ex) {

            ex.printStackTrace();
        }
    }

    StringBuilder contentBuilder = new StringBuilder();

    Features features = null;
    //RepeatFinder finder = null;

//    public TMScoreAnalysis(String pdbCode, String pdbChain) {
//        finder = new RepeatFinder(pdbCode, pdbChain);
//        features = finder.getFeatures();
//    }

    public TMScoreAnalysis(Features features) {
        //this.finder = finder;
        //features = finder.getFeatures();
        this.features = features;
    }

    public void findRepeat() {

        Atom[] atoms = features.getAtoms();
        String strSS = features.getStrSS();

        MutilAlign align = new MutilAlign();
        List<Atom> lAtoms = new ArrayList<Atom>();

        for (Atom atom : atoms)
            lAtoms.add(atom);
        double minRMSD = Double.MAX_VALUE;
        double minGaps = Double.MAX_VALUE;
        double maxScore = Double.MIN_VALUE;

        int maxOptLength = Integer.MIN_VALUE;
        double maxTM = 0.01;//Double.MIN_VALUE;
        double minRMSD_Real = Double.MAX_VALUE;

        /**
         * STARTING
         */
        for (int t = 12; t < 100; t = t + 1) {
//            System.out.println();
            //System.out.println(t);
            /**
             * for each of window size, we take only the minimum RMSD
             */


//        for(int step =0 ; step < t ; step ++) {
//
//            System.out.println();
//            for (int i = step; i < atoms.length - 2 * t; i = i + t) {
//                Atom[] s1 = Fragement.getFragementsofAtoms(atoms, i, i + t - 1);
//                Atom[] s2 = Fragement.getFragementsofAtoms(atoms, i + t, i + 2 * t - 1);
//
//                try {
//                    AFPChain afpChain = pairAlign.pairAlign(s1, s2);
//                    System.out.print(NumberFormatUtils.format(afpChain.getTMScore()) + " ");
//                } catch (StructureException e) {
//                    e.printStackTrace();
//                }
//            }
//        }
//
//
//            break;


            for (int i = 0; i < atoms.length - 2 * t; i = i + 1) {

                try {
                    Atom[] s1 = Fragement.getFragementsofAtoms(atoms, i, i + t - 1);
                    Atom[] s2 = Fragement.getFragementsofAtoms(atoms, i + t, i + 2 * t - 1);

                    String pattern1 = VectorShape.getSSPattern(strSS
                            .substring(i, i + t - 1));
                    String pattern2 = VectorShape.getSSPattern(strSS
                            .substring(i + t, i + 2 * t - 1));
                    if (pattern1.length() >= 2 && pattern2.length() >= 2) {


                        int minL = Math.min(s1.length, s2.length);
                        int maxL = Math.max(s1.length, s2.length);
                        //System.out.println(minL + " and " + maxL);
                        AFPChain afpChain = align.pairAlign(s1, s2);
//                        System.out.print(NumberFormatUtils.format(afpChain.getTMScore()) + ", ");
                        if (afpChain.getOptLength() > 10) {
                            if (pattern1.length() == pattern2.length() && pattern1.length() == 2) {
                                if (pattern1.equals(pattern2)) {
                                    if (maxTM < afpChain.getTMScore()) {
                                        minRMSD = afpChain.getChainRmsd();
                                        maxOptLength = afpChain.getOptLength();
                                        minGaps = afpChain.getGapLen();
                                        maxTM = afpChain.getTMScore();

                                    }
                                }
                            } else {

                                if (maxTM < afpChain.getTMScore()) {
                                    minRMSD = afpChain.getChainRmsd();
                                    maxOptLength = afpChain.getOptLength();
                                    minGaps = afpChain.getGapLen();
                                    maxTM = afpChain.getTMScore();

                                }
                            }


                        }


                    }
                } catch (Exception ex) {
                    ex.printStackTrace();

                }


            }

        }

        /**
         * END LOOP
         */

        /**
         * SYMMETRY
         */



        double symmetryScore = 0.0;


        try {
//            Atom[] ca1 = atoms;
//            Atom[] ca2 =  StructureTools.cloneCAArray(ca1);
//            CeSymm ceSymm = new CeSymm();
//            AFPChain afpChain = ceSymm.pairAlign(ca1, ca2);
//            int symmNr = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
////            System.out.println("Symmetry order of: " + symmNr);
////            System.out.println(afpChain.getTMScore());
//
//            symmetryScore = symmNr + afpChain.getTMScore();
        } catch (Exception e) {
            //e.printStackTrace();
        }
        // write here

//        String content = NumberFormatUtils.format(maxTM)+  "\t" + NumberFormatUtils.format(symmetryScore);
//        if (content.length() > 2 && content.length() < 50) {
//            //System.out.println(content);
//            contentBuilder.append(content);
//            this.tmScore = maxTM;
//            this.ceScore = symmetryScore;
//        }

        this.tmScore = maxTM;
        this.ceScore = symmetryScore;


    }


    public double getTmScore() {
        return tmScore;
    }

    public void setTmScore(double tmScore) {
        this.tmScore = tmScore;
    }

    private double tmScore = 0.0;

    public double getCeScore() {
        return ceScore;
    }

    public void setCeScore(double ceScore) {
        this.ceScore = ceScore;
    }

    private double ceScore = 0.0;
    public String getOutput() {
        return contentBuilder.toString();
    }
}
