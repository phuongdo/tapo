package org.cnrs.crbm.lib.trsfinder;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.math.PeakDetector;
import org.cnrs.crbm.lib.math.Peaks;
import org.cnrs.crbm.lib.raphael.Raphael;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.SuperimposerOutput;
import org.cnrs.crbm.lib.repeats.module.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author pdoviet
 */
public class RMSDSubFinder extends Finder {

    final static boolean debug = false;

    static Logger logger = LoggerFactory.getLogger(RMSDSubFinder.class);

    public void setListWinsize(List<Integer> listWinsize) {
        this.listWinsize = listWinsize;
    }

    List<Integer> listWinsize = new ArrayList<Integer>();

    public void setLstRegions(List<Region> lstRegions) {
        this.lstRegions = lstRegions;
    }

    List<Region> lstRegions = new ArrayList<Region>();

    public RMSDSubFinder(Features features) {
        super(features);
        this.name = "RMSD";

    }

    @Override
    public void findRepeat(Features features) {
        TMEvaluation tmEvaluation = new TMEvaluation();
        double maxTMScore = 0.0;
        for (int noElement : listWinsize) {
            FinderOutput finderOutput = tmEvaluation.findRepeatsBasedOnTMScoreandSecondaryStructure(noElement, lstRegions, features);
            maxTMScore = Math.max(finderOutput.getCombineScore().getTmScore(), maxTMScore);
            if (finderOutput.getRepeats().size() > 0) {
                this.repeats.addAll(finderOutput.getRepeats());
            }
        }
        this.combineScore.setTmScore(maxTMScore);
    }


//    private List<Repeat> extractRepeats(Features features) {
//        TMEvaluation tmEvaluation = new TMEvaluation();
//        List<Repeat> pures = new ArrayList<Repeat>();
//        Atom[] atoms = features.getAtoms();
//
//        double maxTMScore = 0.0;
//        if (atoms.length < 500) {
//            for (int winsize : listWinsize) {
//                FinderOutput finderOutput = tmEvaluation.findRepeatsBasedOnTMScore(winsize, 0, atoms.length - 1, features);
//                maxTMScore = Math.max(finderOutput.getCombineScore().getTmScore(), maxTMScore);
//                pures.addAll(finderOutput.getRepeats());
////                System.out.println(winsize + ":" + finderOutput.getRepeats());
////                System.out.println(maxTMScore);
//            }
//        }
//        // store this maxTM score for SVM
//        this.combineScore.setTmScore(maxTMScore);
//
//        return pures;
//
//    }

    @Deprecated
    private List<Repeat> getPureRMSDFinder(Features features) {
        List<Repeat> pures = new ArrayList<Repeat>();
        Atom[] atoms = features.getAtoms();
        // segmentation
        if (atoms.length < 500) {
            int[] wlist = new int[]{200};
            for (Integer ssize : wlist) {
                for (int i = 0; i < atoms.length; ) {
                    try {
                        int start = i;
                        int end = i + ssize;

                        if (end > atoms.length)
                            end = atoms.length;

                        if ((double) (end - start + 1) / ssize < 0.2)
                            break;

                        for (int winsize : listWinsize) {
                            List<Repeat> lrepeat = AtomFinder.getPureRMSDFinder(features, winsize, start, end);
                            if (lrepeat.size() > 0) pures.addAll(lrepeat);
                        }
                        i = i + ssize;
                    } catch (Exception ex) {
                        ///
                    }

                }
            }
        }

        return pures;

    }

    public static List<Region> generateRegions(String strSS) {
        List<Region> regions = getSSRegion(strSS);
        List<Integer> lstLoopMids = new ArrayList<Integer>();
        lstLoopMids.add(0);
        for (int i = 0; i < regions.size() - 1; i++) {
            Region r1 = regions.get(i);
            Region r2 = regions.get(i + 1);
            int midL = (r1.getEnd() + r2.getStart()) / 2;
            lstLoopMids.add(midL);
        }
        lstLoopMids.add(strSS.length() - 1);
        List<Region> lstRegions = new ArrayList<Region>();
        for (int i = 0; i < lstLoopMids.size() - 1; i++) {
            Region region = new Region(lstLoopMids.get(i), lstLoopMids.get(i + 1) - 1);
            lstRegions.add(region);
        }
        return lstRegions;
    }

    public static void main(String[] args) throws StructureException {
        Superimposer superimposer = new Superimposer();
        String pdbCode = "1koo";
        String pdbChain = "C";
        // generate possible repeat length
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
//        RMSDSubFinder finder = new RMSDSubFinder(repeatFinder.getFeatures());


//        listPosL.add(33);
//        finder.setListWinsize(listPosL);
//        System.out.println("start finding");
//        finder.start();
//        System.out.println("end finding");
        String strSS = repeatFinder.getStrSS();

        System.out.println(strSS);
        List<Region> lstRegions = generateRegions(strSS);
        RMSDSubFinder finder = new RMSDSubFinder(repeatFinder.getFeatures());
        List<Integer> listPosL = new ArrayList<Integer>();
        for (int i = 2; i < lstRegions.size() / 2; i = i + 1) {
            listPosL.add(i);
        }
        finder.setLstRegions(lstRegions);
        finder.setListWinsize(listPosL);
        finder.start();

        for (Repeat repeat : finder.getRepeats()) {
            System.out.println(repeatFinder.showRepeatsByConsole(repeat));
        }

//        int L = 31;
//
//        TMEvaluation tmEvaluation = new TMEvaluation();
//        FinderOutput finderOutput = tmEvaluation.findRepeatsBasedOnTMScore(L, 0, atoms.length, repeatFinder.getFeatures());
//        List<Repeat> repeats = finderOutput.getRepeats();
//        for (Repeat repeat : repeats) {
//            System.out.println(repeatFinder.showRepeatsByConsole(repeat));
//        }


//        StringBuffer buffer = new StringBuffer();
//        for (int i = 0; i < strSS.length(); i++) {
//            if (lstLoopMids.contains(i)) {
//                buffer.append("|");
//            } else buffer.append(".");
//        }
//        System.out.println(buffer.toString());


    }

    public static List<Region> getSSRegion(String secondaryStr) {
        List<Region> regions = new ArrayList<Region>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < secondaryStr.length(); i++) {
            char ch = secondaryStr.charAt(i);
            if (ch == 'B') {
                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'B') {
                    i++;
                }
                i = end = i - 1;
                // save
                regions.add(new Region(start, end));


            } else if (ch == 'H') {

                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'H') {
                    i++;
                }
                i = end = i - 1;

                // check length first

                regions.add(new Region(start, end));

            }


        }

        return regions;
    }

}
