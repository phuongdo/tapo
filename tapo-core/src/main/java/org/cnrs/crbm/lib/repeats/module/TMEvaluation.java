package org.cnrs.crbm.lib.repeats.module;


import msa.Msa;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.conf.ThresholdConfig;
import org.cnrs.crbm.lib.math.Peaks;
import org.cnrs.crbm.lib.multalign.PSim;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.SuperimposerOutput;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by pdoviet on 3/14/2015.
 */
public class TMEvaluation {

    Superimposer superimposer = new Superimposer();
    final static double PATTERN_SECONDARY_THRES = 0.8;

    /**
     * Generator of contiguous segments from L res.
     * <p/>
     * |-----------|------------|------------|
     * L           L              L
     * ---|-----------|------------|---------|
     * L           L            l
     * --------|-----------|------------|----|
     * L            L         ignore too small part
     *
     * @param L
     * @param proteinLength
     * @return
     */
//    public List<Repeat> generateCandidates(int L, int proteinLength) {
//        List<Repeat> candidates = new ArrayList<Repeat>();
//        for (int shift = 0; shift < L; shift = shift + 3) {
//            // extract repeat candidates
//            Repeat candiate = new Repeat();
//            for (int i = 0; i < proteinLength - L; i = i + L) {
//                RepeatContent ru = new RepeatContent(i, i + L - 1);
//                if (ru.size() > L / 2)
//                    candiate.getRepeats().add(ru);
//            }
//            if (candiate.getRepeats().size() > 0)
//                candidates.add(candiate);
//        }
//        return candidates;
//    }

    /**
     * This function aim to generate the TR candidates in a certain region of one protein structure
     * it covert the first function :(
     *
     * @param L
     * @param startRegion
     * @param endRegion
     * @return
     */
    public List<Repeat> generateCandidates(int L, int startRegion, int endRegion) {

        //int proteinLength = endRegion - startRegion + 1;
        List<Repeat> candidates = new ArrayList<Repeat>();
        for (int shift = 0; shift < L; shift = shift + 3) {
            // extract repeat candidates
            Repeat candiate = new Repeat();
            for (int i = startRegion + shift; i < endRegion; i = i + L) {

                int start = i;
                int end = i + L - 1;
                if (end > endRegion - 1)
                    end = endRegion - 1;
                RepeatContent ru = new RepeatContent(start, end);
                if (ru.size() == L)
                    candiate.getRepeats().add(ru);
            }
            if (candiate.getRepeats().size() > 0)
                candidates.add(candiate);
        }
        return candidates;
    }


//    @Deprecated
//    public FinderOutput findRepeatsBasedOnTMScoreV100(int L, int startRegion, int endRegion, Features features) {
//        FinderOutput output = new FinderOutput();
//        List<Repeat> candidates = this.generateCandidates(L, startRegion, endRegion);
//        double maxTMScore = 0.0;
//        for (Repeat candidate : candidates) {
//            try {
//                FinderOutput outCand = this.findRepeatsOfOneSegments(candidate, features);
//                output.getRepeats().addAll(outCand.getRepeats());
//                maxTMScore = Math.max(output.getCombineScore().getTmScore(), outCand.getCombineScore().getTmScore());
//                //System.out.println(maxTMScore);
//
//            } catch (StructureException e) {
//                e.printStackTrace();
//            }
//        }
//
//        output.getCombineScore().setTmScore(maxTMScore);
//
//        return output;
//    }


    public FinderOutput findRepeatsBasedOnTMScoreandSecondaryStructure(int L, List<Region> regions, Features features) {
        FinderOutput output = new FinderOutput();
        Atom[] atoms = features.getAtoms();
        String strSS = features.getStrSS();
        Atom[] atomSetTemplate = null;
        Atom[] atomSetQuery = null;
        double maxTMScore = 0.0;
        int startRegion = 0;
        int endRegion = regions.size() - 1;
        List<Double> datas = new ArrayList<Double>();
        for (int i = startRegion; i < endRegion - 2 * L; i++) {
            atomSetTemplate = Fragement.getFragementsofAtoms(atoms, regions.get(i).getStart(), regions.get(i + L - 1).getEnd());
            String pattern1 = VectorShape.getSSPattern(strSS
                    .substring(regions.get(i).getStart(), regions.get(i + L - 1).getEnd()));
            atomSetQuery = Fragement.getFragementsofAtoms(atoms, regions.get(i + L).getStart(), regions.get(i + L + L - 1).getEnd());
            String pattern2 = VectorShape.getSSPattern(strSS
                    .substring(regions.get(i + L).getStart(), regions.get(i + L + L - 1).getEnd()));
            boolean check = (Math.abs(atomSetTemplate.length - atomSetQuery.length) / ((atomSetQuery.length + atomSetTemplate.length) / 2)) < 0.2;
            /**
             * from version 1.1.3-unstable (coming soon)
             */
            String[] copies = new String[]{pattern1, pattern2};
            double scoreSS = this.scoreSSPattern(copies);
            if (pattern1.length() >= 2 && pattern2.length() >= 2 && check && scoreSS >= PATTERN_SECONDARY_THRES) {
                try {
                    SuperimposerOutput superOuput = superimposer.compareTwoStructures(atomSetTemplate, atomSetQuery);
                    datas.add(superOuput.getTmScore());
                    maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
                } catch (Exception ex) {
//                    ex.printStackTrace();
                    datas.add(0.0);
                }

            } else {
                datas.add(0.0);
            }
        }

        double[] signals = new double[datas.size()];
        for (int i = 0; i < signals.length; i++) {
            signals[i] = datas.get(i);
        }
        // finding the best
        LinkedList<Integer> peaks = Peaks.findPeaks(signals, L, ThresholdConfig.TM_THRES);
        //
        //System.out.println(peaks);
        for (int peak : peaks) {
            int i = peak + startRegion;
            Repeat aRepeat = new Repeat();
            atomSetTemplate = Fragement.getFragementsofAtoms(atoms, regions.get(i).getStart(), regions.get(i + L - 1).getEnd());
            aRepeat.getRepeats().add(new RepeatContent(regions.get(i).getStart(), regions.get(i + L - 1).getEnd()));
            String pattternTemplate = VectorShape.getSSPattern(strSS
                    .substring(regions.get(i).getStart(), regions.get(i + L - 1).getEnd()));
            /**
             * scan in the right
             * ----------xxxxxx[----------]----
             */
            //scan right;
            for (int j = i + L; j < endRegion - L; j = j + L) {
                atomSetQuery = Fragement.getFragementsofAtoms(atoms, regions.get(j).getStart(), regions.get(j + L - 1).getEnd());
                // compare
                String patternQuery = VectorShape.getSSPattern(strSS
                        .substring(regions.get(j).getStart(), regions.get(j + L - 1).getEnd()));
                if (pattternTemplate.length() >= 2 && pattternTemplate.equals(patternQuery)) {
                    try {
                        SuperimposerOutput superOuput = superimposer.compareTwoStructures(atomSetTemplate, atomSetQuery);
                        maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
                        if (superOuput.isTRs()) {
                            aRepeat.getRepeats().add(
                                    new RepeatContent(regions.get(j).getStart(), regions.get(j + L - 1).getEnd()));
                        } else {
                            break;
                        }

                    } catch (Exception ex) {
                        ex.printStackTrace();
                        break;
                    }
                } else {
                    //fuck
                    break;
                }

            }// end scan right

            /**
             * scan in the left
             * ----------[----------]xxxxxx----
             */

            for (int j = i - L; j >= startRegion; j = j - L) {
                atomSetQuery = Fragement.getFragementsofAtoms(atoms, regions.get(j).getStart(), regions.get(j + L - 1).getEnd());
                // compare
                String patternQuery = VectorShape.getSSPattern(strSS
                        .substring(regions.get(j).getStart(), regions.get(j + L - 1).getEnd()));
                if (pattternTemplate.length() >= 2 && pattternTemplate.equals(patternQuery)) {
                    try {
                        SuperimposerOutput superOuput = superimposer.compareTwoStructures(atomSetTemplate, atomSetQuery);
                        maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
                        if (superOuput.isTRs()) {
                            aRepeat.getRepeats().add(
                                    new RepeatContent(regions.get(j).getStart(), regions.get(j + L - 1).getEnd()));
                        } else {
                            break;
                        }

                    } catch (Exception ex) {
                        ex.printStackTrace();
                        break;
                    }
                } else {
                    //fuck
                    break;
                }

            }// end scan left

            if (aRepeat.getRepeats().size() >= 2) {
                aRepeat.sortByPosition();
                output.getRepeats().add(aRepeat);
            }
        }

        output.getCombineScore().setTmScore(maxTMScore);
        return output;
    }


    public FinderOutput findRepeatsBasedOnTMScore(int L, int startRegion, int endRegion, Features features) {
        FinderOutput output = new FinderOutput();
        Atom[] atoms = features.getAtoms();
        String strSS = features.getStrSS();
        Atom[] atomSetTemplate = null;
        Atom[] atomSetQuery = null;
        double maxTMScore = 0.0;
        List<Double> datas = new ArrayList<Double>();
        for (int i = startRegion; i <= endRegion - 2 * L + 1; i++) {
            atomSetTemplate = Fragement.getFragementsofAtoms(atoms, i, i + L - 1);
            String pattern1 = VectorShape.getSSPattern(strSS
                    .substring(i, i + L - 1));
            atomSetQuery = Fragement.getFragementsofAtoms(atoms, i + L, i + L + L - 1);
            String pattern2 = VectorShape.getSSPattern(strSS
                    .substring(i + L, i + L + L - 1));
            /**
             * from version 1.1.3-unstable (coming soon)
             */

            String[] copies = new String[]{pattern1, pattern2};
            double scoreSS = this.scoreSSPattern(copies);
            if (pattern1.length() >= 2 && pattern2.length() >= 2 && scoreSS >= PATTERN_SECONDARY_THRES) {
                try {
                    SuperimposerOutput superOuput = superimposer.compareTwoStructures(atomSetTemplate, atomSetQuery);
                    datas.add(superOuput.getTmScore());
                    maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
                } catch (Exception ex) {
                    datas.add(0.0);
                }

            } else {
                datas.add(0.0);
            }
        }

        double[] signals = new double[datas.size()];
        for (int i = 0; i < signals.length; i++) {
            signals[i] = datas.get(i);
        }


        // finding the best
        LinkedList<Integer> peaks = Peaks.findPeaks(signals, L, ThresholdConfig.TM_THRES);
//        System.out.println(peaks);

        for (int peak : peaks) {
            int i = peak + startRegion;
            Repeat aRepeat = new Repeat();
            atomSetTemplate = Fragement.getFragementsofAtoms(atoms, i, i + L - 1);
            aRepeat.getRepeats().add(new RepeatContent(i, i + L - 1));
            String pattternTemplate = VectorShape.getSSPattern(strSS
                    .substring(i, i + L - 1));
            /**
             * scan in the right
             * ----------xxxxxx[----------]----
             */
            //scan right;
            for (int j = i + L; j <= endRegion - L + 1; j = j + L) {
                atomSetQuery = Fragement.getFragementsofAtoms(atoms, j, j + L - 1);
                // compare
                String patternQuery = VectorShape.getSSPattern(strSS
                        .substring(j, j + L - 1));
                if (pattternTemplate.length() >= 2 && pattternTemplate.equals(patternQuery)) {
                    try {
                        SuperimposerOutput superOuput = superimposer.compareTwoStructures(atomSetTemplate, atomSetQuery);
                        maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
                        if (superOuput.isTRs()) {
                            aRepeat.getRepeats().add(
                                    new RepeatContent(j, j + L - 1));
                        } else {
                            break;
                        }

                    } catch (Exception ex) {
                        ex.printStackTrace();
                        break;
                    }
                } else {
                    //fuck
                    break;
                }

            }// end scan right

            /**
             * scan in the left
             * ----------[----------]xxxxxx----
             */

            for (int j = i - L; j >= startRegion; j = j - L) {
                atomSetQuery = Fragement.getFragementsofAtoms(atoms, j, j + L - 1);
                // compare
                String patternQuery = VectorShape.getSSPattern(strSS
                        .substring(j, j + L - 1));
                if (pattternTemplate.length() >= 2 && pattternTemplate.equals(patternQuery)) {
                    try {
                        SuperimposerOutput superOuput = superimposer.compareTwoStructures(atomSetTemplate, atomSetQuery);
                        maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
                        if (superOuput.isTRs()) {
                            aRepeat.getRepeats().add(
                                    new RepeatContent(j, j + L - 1));
                        } else {
                            break;
                        }

                    } catch (Exception ex) {
                        ex.printStackTrace();
                        break;
                    }
                } else {
                    //fuck
                    break;
                }

            }// end scan left

            if (aRepeat.getRepeats().size() >= 2) {
                aRepeat.sortByPosition();
                output.getRepeats().add(aRepeat);
            }
        }

        output.getCombineScore().setTmScore(maxTMScore);
        return output;
    }

    public FinderOutput findRepeatsOfOneSegments(Repeat queryRepeat, Features features) throws StructureException {
        FinderOutput output = new FinderOutput();
        List<Repeat> repeats = new ArrayList<Repeat>();
        Atom[] atoms = features.getAtoms();
        String strSS = features.getStrSS();
        double maxTMScore = 0.0;
        List<RepeatContent> candidates = queryRepeat.getRepeats();
        int[] spaces = new int[candidates.size() - 1];

        for (int i = 0; i < candidates.size() - 1; i++) {
            Atom[] query = Fragement.getFragementsofAtoms(atoms, candidates.get(i).getStart(), candidates.get(i).getEnd());
            Atom[] target = Fragement.getFragementsofAtoms(atoms, candidates.get(i + 1).getStart(), candidates.get(i + 1).getEnd());
            String pattern1 = VectorShape.getSSPattern(strSS
                    .substring(candidates.get(i).getStart(), candidates.get(i).getEnd()));
            String pattern2 = VectorShape.getSSPattern(strSS
                    .substring(candidates.get(i + 1).getStart(), candidates.get(i + 1).getEnd()));

//            int start1 = candidates.get(i).getStart();
//            int end1 = candidates.get(i).getEnd();
//            int start2 = candidates.get(i + 1).getStart();
//            int end2 = candidates.get(i + 1).getEnd();
//            boolean checkLoop = strSS.charAt(candidates.get(i).getStart()) == '-' &&
//                    strSS.charAt(candidates.get(i).getEnd()) == '-' && strSS.charAt(candidates.get(i + 1).getEnd()) == '-';
            //if (checkLoop && ((pattern1.equals(pattern2) && pattern1.length() == 2) || (pattern1.length() >= 2 && pattern2.length() >= 2))) {
            if (pattern1.equals(pattern2) && pattern1.length() >= 2) {
                SuperimposerOutput superOuput = superimposer.compareTwoStructures(query, target);
                maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
                if (superOuput.isTRs())
                    spaces[i] = 1;
            }
        }
        List<Integer> breakPoints = new ArrayList<Integer>();
        breakPoints.add(0);
        for (int i = 0; i < spaces.length; i++) {
            if (spaces[i] == 1) {
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

        for (RepeatContent rc : candidates.subList(start, end + 1)) {
            splitRepeat.getRepeats().add(rc);
        }
        // splitRepeat.setRepeats(contents.subList(start, end + 1));

        if (splitRepeat.getRepeats().size() >= 2)
            repeats.add(splitRepeat);

        for (int i = 1; i < breakPoints.size() - 1; i++) {

            // System.out.println(breakPoints.get(i) + 1 + "-"
            // + breakPoints.get(i + 1));

            start = breakPoints.get(i) + 1;
            end = breakPoints.get(i + 1);

            splitRepeat = new Repeat();
            // splitRepeat.setRepeats(contents.subList(start, end + 1));
            for (RepeatContent rc : candidates.subList(start, end + 1)) {
                splitRepeat.getRepeats().add(rc);
            }

            if (splitRepeat.getRepeats().size() >= 2)
                repeats.add(splitRepeat);

        }

        output.setRepeats(repeats);
        output.getCombineScore().setTmScore(maxTMScore);
        return output;
    }


    public double getMaxTMScore(Repeat queryRepeat, Features features) throws StructureException {
        Atom[] atoms = features.getAtoms();
        String strSS = features.getStrSS();
        double maxTMScore = 0.0;
        List<RepeatContent> candidates = queryRepeat.getRepeats();
        for (int i = 0; i < candidates.size() - 1; i++) {
            Atom[] query = Fragement.getFragementsofAtoms(atoms, candidates.get(i).getStart(), candidates.get(i).getEnd());
            Atom[] target = Fragement.getFragementsofAtoms(atoms, candidates.get(i + 1).getStart(), candidates.get(i + 1).getEnd());
            String pattern1 = VectorShape.getSSPattern(strSS
                    .substring(candidates.get(i).getStart(), candidates.get(i).getEnd()));
            String pattern2 = VectorShape.getSSPattern(strSS
                    .substring(candidates.get(i + 1).getStart(), candidates.get(i + 1).getEnd()));
            if (pattern1.equals(pattern2) && pattern1.length() >= 2) {
                SuperimposerOutput superOuput = superimposer.compareTwoStructures(query, target);
                maxTMScore = Math.max(maxTMScore, superOuput.getTmScore());
            }
        }

        return maxTMScore;
    }


    /**
     * An improvement of AtomFinder class. This function combines 2 advantages of Extension Module and AtomFinderModule
     *
     * @param queryRepeat
     * @param features
     * @return
     * @throws StructureException
     */
    @Deprecated
    public List<Repeat> tmEvaluationRepeat(Repeat queryRepeat, Features features) throws StructureException {

        List<Repeat> repeats = new ArrayList<Repeat>();
        Atom[] atoms = features.getAtoms();
        List<RepeatContent> candidates = queryRepeat.getRepeats();
        int[] spaces = new int[candidates.size() - 1];
        for (int i = 0; i < candidates.size() - 1; i++) {

            Atom[] query = Fragement.getFragementsofAtoms(atoms, candidates.get(i).getStart(), candidates.get(i).getEnd());
            Atom[] target = Fragement.getFragementsofAtoms(atoms, candidates.get(i + 1).getStart(), candidates.get(i + 1).getEnd());
            if (superimposer.isTheSameStructure(query, target))
                spaces[i] = 1;
            else spaces[i] = 0;
        }

        List<Integer> breakPoints = new ArrayList<Integer>();
        breakPoints.add(0);
        for (int i = 0; i < spaces.length; i++) {
            if (spaces[i] == 1) {
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

        for (RepeatContent rc : candidates.subList(start, end + 1)) {
            splitRepeat.getRepeats().add(rc);
        }
        // splitRepeat.setRepeats(contents.subList(start, end + 1));

        if (splitRepeat.getRepeats().size() >= 2)
            repeats.add(splitRepeat);

        for (int i = 1; i < breakPoints.size() - 1; i++) {

            // System.out.println(breakPoints.get(i) + 1 + "-"
            // + breakPoints.get(i + 1));

            start = breakPoints.get(i) + 1;
            end = breakPoints.get(i + 1);

            splitRepeat = new Repeat();
            // splitRepeat.setRepeats(contents.subList(start, end + 1));
            for (RepeatContent rc : candidates.subList(start, end + 1)) {
                splitRepeat.getRepeats().add(rc);
            }

            if (splitRepeat.getRepeats().size() >= 2)
                repeats.add(splitRepeat);

        }
        return repeats;
    }


    public double scoreSSPattern(String[] copies) {

        //String[] copies = new String[]{"HB", "BBH", "HB"};
//        String[] copies = new String[]{"HB", "BBH", "HB"};
        double score = 0.0;
        try {
            Msa ms = new Msa(copies, "");
            LinkedList<String> aligned = ms.buildAlignment();
            List<String> msa = new ArrayList<String>();
            for (String a : aligned) {
                msa.add(a);
            }
            PSim pSim = new PSim(msa);
            score = pSim.getSimilarity();
        } catch (Exception ex) {
        }
        return score;


    }

}
