package org.cnrs.crbm.lib.trsfinder;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.conf.ConfigUtil;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.CombineScore;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatExtend;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

public abstract class Finder {

    static Logger logger = LoggerFactory.getLogger(Finder.class);
    Superimposer superimposer = new Superimposer();
    TMEvaluation tmEvaluation = new TMEvaluation();

    protected Features features = null;
    protected List<Repeat> repeats = new ArrayList<Repeat>();
    protected String name = "Finder defaults";
    protected boolean isTRs = false;
    protected int nrOfTRs = 0;// new feature
    protected Atom[] atoms = null;
    RepeatExtend repeatExtend = null;
    protected double theBestScore = 0.0;
    protected CombineScore combineScore = new CombineScore();

    public boolean isTRs() {

        return this.isTRs;
    }

    public int getNrOfTRs() {

        if (repeats.size() > 0) {
            int tmpNrOfTRs = 0;
            for (Repeat repeat : repeats) {
                tmpNrOfTRs = repeat.getRepeats().size();
                this.nrOfTRs = Math.max(tmpNrOfTRs, nrOfTRs);

            }

        }

        return this.nrOfTRs;
    }

    protected List<Repeat> splitRepeat(Repeat repeat) {
        List<Repeat> repeats_copy = new ArrayList<Repeat>();
        double avgLengh = repeat.getAvgLength();
        List<RepeatContent> contents = repeat.getRepeats();
        int[] spaces = new int[contents.size() - 1];

        for (int i = 0; i < contents.size() - 1; i++) {
            spaces[i] = contents.get(i + 1).getStart()
                    - contents.get(i).getEnd();
        }

        List<Integer> breakPoints = new ArrayList<Integer>();
        breakPoints.add(0);
        for (int i = 0; i < spaces.length; i++) {

            if (spaces[i] < (double) 0.6 * avgLengh) {

            } else {
                breakPoints.add(i);
                // System.out.println(i);
            }

        }

        breakPoints.add(spaces.length);

        // first
        // System.out.println(breakPoints.get(0) + "-" + breakPoints.get(1));

        // REMOVE the repeat has only one segment
        int start = breakPoints.get(0);
        int end = breakPoints.get(1);

        Repeat splitRepeat = new Repeat();

        for (RepeatContent rc : contents.subList(start, end + 1)) {
            splitRepeat.getRepeats().add(rc);
        }
        // splitRepeat.setRepeats(contents.subList(start, end + 1));

        if (splitRepeat.getRepeats().size() >= 2)
            repeats_copy.add(splitRepeat);

        for (int i = 1; i < breakPoints.size() - 1; i++) {

            // System.out.println(breakPoints.get(i) + 1 + "-"
            // + breakPoints.get(i + 1));

            start = breakPoints.get(i) + 1;
            end = breakPoints.get(i + 1);

            splitRepeat = new Repeat();
            // splitRepeat.setRepeats(contents.subList(start, end + 1));
            for (RepeatContent rc : contents.subList(start, end + 1)) {
                splitRepeat.getRepeats().add(rc);
            }

            if (splitRepeat.getRepeats().size() >= 2)
                repeats_copy.add(splitRepeat);

        }

        // last

        // if (split_repeat.getRepeats().size() >= 2)
        // repeats_copy.add(split_repeat);

        return repeats_copy;

    }

    // public static void main(String[] args) {
    // Repeat repeat = new Repeat();
    // repeat.getRepeats().add(new RepeatContent(1, 20));
    // repeat.getRepeats().add(new RepeatContent(21, 40));
    // repeat.getRepeats().add(new RepeatContent(41, 60));
    // repeat.getRepeats().add(new RepeatContent(90, 110));
    // repeat.getRepeats().add(new RepeatContent(111, 130));
    //
    // repeat.getRepeats().add(new RepeatContent(596, 619));
    // repeat.getRepeats().add(new RepeatContent(621, 639));
    //
    // List<Repeat> repeats = splitRepeat(repeat);
    // for (Repeat arepeat : repeats) {
    // System.out.println("==");
    // System.out.println(arepeat);
    // }
    // }

    private boolean checkTRs() {

        /**
         * AFTER 26 May 2015
         * we remove this part out of taking consideration
         */

//		List<Repeat> repeats_copy = new ArrayList<Repeat>();
//		if (repeats.size() > 0) {
//			for (Repeat repeat : repeats) {
//				// double avgLengh = repeat.getAvgLength();
//				List<Repeat> splits = this.splitRepeat(repeat);
//				int tmpNrOfTRs = 0;
//				// int numberOfTRs = 0;
//				for (Repeat split : splits) {
//					tmpNrOfTRs = split.getRepeats().size();
//					this.nrOfTRs = Math.max(tmpNrOfTRs, nrOfTRs);
//					repeats_copy.add(split);
//					// System.out.println(split);
//				}
//			}
//
//		}
//		this.repeats = repeats_copy;
        if (repeats.size() > 0)
            return true;
        else
            return false;

    }

    private void filterTRs() {

        List<Repeat> repeats_copy = new ArrayList<Repeat>();
        String strSS = features.getStrSS();
        // int proteinSize = strSS.length();
        if (repeats.size() > 0) {

            Atom[] atomSet1 = null;
            Atom[] atomSet2 = null;
            for (Repeat repeat : repeats) {


                List<RepeatContent> contents = repeat.getRepeats();
                int[] labels = new int[contents.size()];
                for (int k = 0; k < contents.size(); k++) {
                    int start1 = contents.get(k).getStart();
                    int end1 = contents.get(k).getEnd();

                    int match = 0;
                    atomSet1 = Fragement.getFragementsofAtoms(
                            features.getAtoms(), start1, end1);
                    for (int i = 0; i < contents.size(); i++) {
                        int start2 = contents.get(i).getStart();
                        int end2 = contents.get(i).getEnd();
                        try {

                            atomSet2 = Fragement.getFragementsofAtoms(
                                    features.getAtoms(), start2, end2);

                            String pattern1 = VectorShape.getSSPattern(strSS
                                    .substring(start1, end1 + 1));
                            String pattern2 = VectorShape.getSSPattern(strSS
                                    .substring(start2, end2 + 1));

//                            String annotation = ProteinAnnotate.anotation(
//                                    features, start2, end2);


                            // double leng = repeat.getAvgLength();
                            // System.out.println(leng + ":" + proteinSize);
                            // if (atomSet1.length < 40
                            // && leng < (double) 0.33 * proteinSize)
                            // check = pattern1.equals(pattern2)
                            // && pattern1.length() >= 2;
                            boolean check = true;
                            if (atomSet1.length < 25) {
                                check = pattern1.equals(pattern2) && (pattern1.length() > 1 || pattern1.length() == 0);
                            } else if (atomSet1.length < 40) {
                                if (pattern1.length() > 1 && pattern1.length() <= 3)
                                    check = pattern1.equals(pattern2);
                                else if (pattern1.length() == 1 && pattern1.equals(pattern2))
                                    check = false;
                                else
                                    check = true;
                            }
//                            if (pattern1.length() == 1 && pattern2.length() == 1)
//                                check = false;
                            if (check
                                //  && !annotation.contains("MA")
                                //&& !annotation.contains("LH")
//                                    && superimposer.isTheSameContactMap(
//                                    atomSet1, start1, atomSet2, start2,
//                                    features.getCmHistos())
                                    ) {
                                match++;
                            }
                        } catch (Exception e) {
                            // TODO Auto-generated catch block
                            // System.out.println(start2 + "-" + end2);
                            logger.error(e.getMessage() + " with pdb "
                                    + features.getPdbCode() + "_"
                                    + features.getPdbChain());
                        }

                    }

                    if ((double) match / contents.size() > 0.5)
                        labels[k] = 1;

                }

                Repeat repeat_copy = new Repeat();
                for (int i = 0; i < contents.size(); i++) {

                    if (labels[i] == 1)
                        repeat_copy.getRepeats().add(contents.get(i));
                }

                if (repeat_copy.getRepeats().size() >= 2)
                    repeats_copy.add(repeat_copy);

            }

        }
        this.repeats = repeats_copy;

    }

    public Finder(Features features) {
        this.features = features;
        atoms = features.getAtoms();

    }

    /**
     * Starting process. This process of finding tandem repeats only starting when
     * this method is called.
     */

    public void start() {
        try {
//            System.out.println(this.getName() + " is starting...");
            this.findRepeat(features);
            String calQA = ConfigUtil.getInstance().getProperty("CALS_QA");
            if (calQA.equals("on")) {
                repeatExtend = new RepeatExtend(features);
//            this.filterTRs();
                this.isTRs = this.checkTRs();
                this.extendingTRs(this);
                // remove duplicate
                this.removeDuplicate();
                // assign QA
                this.assignQA();
//                System.out.println(this.getName() + " is done!");

                /**
                 * updated from version 1.1.0
                 */
                this.calucateMaxTMScore();
            }
        } catch (Exception ex) {
            // do nothing
//            ex.printStackTrace();
            logger.error(ex.getMessage() + "Finder() with pdb "
                    + features.getPdbCode() + "_" + features.getPdbChain());
        }
    }


    private void calucateMaxTMScore() {

        double maxTMScore = 0.0;
        for (Repeat repeat : repeats) {
            try {
                maxTMScore = Math.max(maxTMScore, tmEvaluation.getMaxTMScore(repeat, features));
            } catch (Exception ex) {

            }
        }
        this.combineScore.setTmScore(Math.max(maxTMScore, combineScore.getTmScore()));
    }

    private void removeDuplicate() {

        List<Repeat> newListRepeat = new ArrayList<Repeat>();
        for (Repeat repeat : repeats) {
            if (!newListRepeat.contains(repeat)) {
                newListRepeat.add(repeat);
            }
        }
        this.repeats = newListRepeat;

    }

    private void assignQA() {
        for (Repeat repeat : repeats) {
            StringBuffer buffer = new StringBuffer();
            for (RepeatContent content : repeat.getRepeats()) {
                buffer.append(content.getStart() + "-" + content.getEnd() + ";");
            }
            String align = buffer.toString();
            align = align.substring(0, align.length() - 1);
            MutilAlign mutilAlign = new MutilAlign(features, align);
            repeat.setScore(mutilAlign.getQAScore());
//            double a = 1;
//            double b = 1;
//            repeat.setScore1(mutilAlign.getScore1());
//            repeat.setScore2(mutilAlign.getScore2());
//            repeat.setScore3(mutilAlign.getScore3());
//            repeat.setScore4(mutilAlign.getScore4());

            //repeat.setRankScore(repeat.getScore() * Math.pow(1 - 1 / (1 + repeat.getRepeats().size()), a) * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));
//            repeat.setRankScore(repeat.getScore() * repeat.getRepeats().size() * Math.pow(1 - 1 / (1 + repeat.getAvgLength()), b));
            //repeat.setRankScore(mutilAlign.getQAScore());

//            System.out.println(repeat.getAvgLength() + "\t" + repeat.getScore());
        }


    }


    public void extendingTRs(Finder finder) {

        List<Repeat> repeats = finder.getRepeats();
        // List<Repeat> repeats_copy = new ArrayList<Repeat>();
        for (Repeat repeat : repeats) {
            try {
                List<RepeatContent> templateTRs = new ArrayList<RepeatContent>();
                for (RepeatContent r : repeat.getRepeats()) {
                    templateTRs.add(r);
                }

                repeatExtend.extendRepeatRight(repeat, templateTRs);
                repeatExtend.extendRepeatLeft(repeat, templateTRs);

            } catch (Exception ex) {
//                    System.out.println(this.getName());
//                    ex.printStackTrace();
            }

            // repeats_copy.add(repeat);

        }


    }


    /**
     * 22/05 Must be fixed the repeats containing only one secondary structure
     *
     * @param repeat
     * @return
     */
    public boolean isSimilaryStructuralRepeat(Repeat repeat) {

        boolean istandem = false;

        String strSS = features.getStrSS();

        List<RepeatContent> contents = repeat.getRepeats();
        for (int i = 0; i < contents.size() - 1; i++) {
            int start1 = contents.get(i).getStart();
            int end1 = contents.get(i).getEnd();
            int start2 = contents.get(i + 1).getStart();
            int end2 = contents.get(i + 1).getEnd();
            try {

                // pattern checking
//				String pattern1 = VectorShape.getSSPattern(strSS.substring(
//						start1, end1 + 1));
//				String pattern2 = VectorShape.getSSPattern(strSS.substring(
//						start2, end2 + 1));

//				if (pattern1.equals(pattern2)) {

                Atom[] atomSet1 = Fragement.getFragementsofAtoms(
                        features.getAtoms(), start1, end1);
                Atom[] atomSet2 = Fragement.getFragementsofAtoms(
                        features.getAtoms(), start2, end2);

                if (superimposer.isTheSameStructure(atomSet1, atomSet2)) {
                    istandem = true;
                    break;
                }
//				}
            } catch (StructureException e) {
                // TODO Auto-generated catch block
                System.out.println(start2 + "-" + end2);
                // e.printStackTrace();
            }

        }

        return istandem;

    }

    public abstract void findRepeat(Features features);

    public List<Repeat> getRepeats() {
        return repeats;
    }

    public void setRepeats(List<Repeat> repeats) {
        this.repeats = repeats;
    }

    public String getName() {
        return this.name;
    }

    public CombineScore getCombineScore() {
        return combineScore;
    }

    public void setCombineScore(CombineScore combineScore) {
        this.combineScore = combineScore;
    }


}
