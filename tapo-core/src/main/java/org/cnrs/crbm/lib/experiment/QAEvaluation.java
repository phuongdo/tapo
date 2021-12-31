package org.cnrs.crbm.lib.experiment;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.*;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 1/8/2015.
 * We want to test how many times will be spent. by using only this
 * RMSD method. See Algorithm 1 in our report.
 */
public class QAEvaluation {

    public static void main(String[] args) {
        String pdbCode = "1fnx";
        String pdbChain = "H";
        if (args.length > 0) {
            pdbCode = args[0];
            pdbChain = args[1];
        }
        try {
            QAEvaluation evaluation = new QAEvaluation(pdbCode, pdbChain);
            evaluation.findRepeat();
            System.out.println(evaluation.getOutput());
            //System.exit(0);

        } catch (Exception ex) {

            ex.printStackTrace();
        }
    }

    StringBuilder contentBuilder = new StringBuilder();

    Features features = null;
    RepeatFinder finder = null;

    public QAEvaluation(String pdbCode, String pdbChain) {
        finder = new RepeatFinder(pdbCode, pdbChain);
        features = finder.getFeatures();
    }

    public void findRepeat() {
        List<Repeat> allTRs = new ArrayList<Repeat>();

        //generate TRs
        Atom[] atoms = features.getAtoms();
        String strSS = features.getStrSS();
        String pdbCode = features.getPdbCode();
        String pdbChain = features.getPdbChain();
        for (int nTRs = 2; nTRs <= 5; nTRs++) {
            for (int winsize = 10; winsize < 90; winsize = winsize + 3) {
                for (int i = 0; i < atoms.length - (winsize * nTRs); i = i + 5) {
                    Repeat repeat = new Repeat();
                    // build a segment with at least 2 secondary structure.
                    for (int n = 0; n < nTRs; n++) {
                        int start = i + n * winsize;
                        int end = i + (n + 1) * winsize - 1;
                        String pattern1 = VectorShape.getSSPattern(strSS
                                .substring(start, end));
                        if (pattern1.length() >= 2)
                            repeat.getRepeats().add(new RepeatContent(start, end));
                    }
                    if (repeat.getRepeats().size() >= 2)
                        allTRs.add(repeat);

                }
            }
        }


        // ranking
        try {

            double bestQA = -1;
            double bestQAOld = -1;
            double bestScore1 = -1;
            double bestScore2 = -1;
            double bestScore3 = -1;
            double bestScore4 = -1;

            double bestRank = -1;
            double bestRank2 = -1;

            for (Repeat repeat : allTRs) {
                StringBuffer buffer = new StringBuffer();
                for (RepeatContent content : repeat.getRepeats()) {
                    buffer.append(content.getStart() + "-" + content.getEnd() + ";");
                }

                String align = buffer.toString();
                align = align.substring(0, align.length() - 1);

                MutilAlign mutilAlign = new MutilAlign(finder.getFeatures(), align);
                repeat.setScore(mutilAlign.getQAScore());
                repeat.setScore1(mutilAlign.getScore1());
                repeat.setScore2(mutilAlign.getScore2());

                double b = 1;
                repeat.setRankScore(repeat.getScore() * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));

                if (bestQA < repeat.getScore()) {
                    bestQA = Math.max(repeat.getScore(), bestQA);
                    bestScore1 = repeat.getScore1();
                    bestScore2 = repeat.getScore2();
                    bestScore3 = mutilAlign.getScore3();
                    bestScore4 = mutilAlign.getScore4();
                    bestRank = repeat.getRankScore();
                    bestQAOld = ((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3);
                    bestRank2 = ((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3) * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b);
//
//                    bestScore1 = Math.max(repeat.getScore1(), bestScore1);
//                    bestScore2 = Math.max(repeat.getScore2(), bestScore2);
//                    bestScore3 = Math.max(mutilAlign.getScore3(), bestScore3);
//                    bestScore4 = Math.max(mutilAlign.getScore4(), bestScore4);
//                    bestRank = Math.max(repeat.getRankScore(), bestRank);
//                    bestQAOld = Math.max(((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3), bestQAOld);
//                    bestRank2 = Math.max(((mutilAlign.getScore4() + (mutilAlign.getScore3() + repeat.getScore2()) / 2) / 3) * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b), bestRank2)
                    ;


                    //repeat.setRankScore(repeat.getScore() * Math.pow(1 - 1 / (1 + repeat.getRepeats().size()), a) * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));
                    //System.out.println(bestScore2);
                }
            }

            contentBuilder.append(pdbCode + "_" + pdbChain + "\t" + NumberFormatUtils.format(bestQAOld) + "\t" + NumberFormatUtils.format(bestQA) + "\t" + NumberFormatUtils.format(bestScore1) + "\t" + NumberFormatUtils.format(bestScore2) + "\t" + NumberFormatUtils.format(bestScore3) + "\t" + NumberFormatUtils.format(bestScore4) + "\t" + NumberFormatUtils.format(bestRank) + "\t" + NumberFormatUtils.format(bestRank2));
        } catch (Exception ex) {

        }


    }

    public String getOutput() {
        return contentBuilder.toString();
    }

}
