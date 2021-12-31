package org.cnrs.crbm.lib.repeats;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.conf.ThresholdConfig;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class RepeatExtend {

    public static void main(String[] args) throws StructureException {
//        RepeatFinder finder = new RepeatFinder("1lxa", "A");
//        Atom[] atoms = finder.getAtoms();
//        int start = 43;
//        int l = 18;

        RepeatFinder finder = new RepeatFinder("1dce", "A");
        Atom[] atoms = finder.getAtoms();
        int start = 124;
        int l = 36;

        Atom[] atomTemp = Fragement.getFragementsofAtoms(atoms,
                start, start + l - 1);
        MutilAlign mutilAlign = new MutilAlign();


        for (int w = l; w < l * 2; w++) {

            Atom[] atomCompare = Fragement.getFragementsofAtoms(atoms,
                    start + l, start + l + w);
            AFPChain afpChain = mutilAlign.pairAlign(atomTemp, atomCompare);

            System.out.println(start + l + w + "\t" + NumberFormatUtils.format(afpChain.getTMScore()));
        }


    }

    MutilAlign mutilAlign = new MutilAlign();
    //Superimposer superimposer = new Superimposer();
    private Atom[] atoms = null;
    String strSS;
    //static final double TM_SCORE_THRESH = 0.4;
    double[] cmhistos = null;

    public RepeatExtend(Features features) {

        this.atoms = features.getAtoms();
        this.strSS = features.getStrSS();
        this.cmhistos = features.getCmHistos();
    }

    public void extendRepeatRight(Repeat repeat, List<RepeatContent> templateTRs) {

        // doing something here
        // int start = repeat.getStart();
        int end = repeat.getEnd();
        // System.out.println(end);
        int trsLen = (int) repeat.getAvgLength();
        // take the last unit as template
        RepeatContent lastUnit = repeat.getRepeats().get(repeat.getRepeats().size() - 1);
        Atom[] atomTemp = Fragement.getFragementsofAtoms(atoms,
                lastUnit.getStart(), lastUnit.getEnd());

        List<Double> points = new ArrayList<Double>();
        for (int window = trsLen; window < atoms.length - end
                && window < trsLen * 2; window++) {
            try {
                Atom[] atomCompare = Fragement.getFragementsofAtoms(atoms,
                        end + 1, end + window);

                AFPChain afpChain = mutilAlign.pairAlign(atomTemp, atomCompare);
                double tmScore = afpChain.getTMScore();
                points.add(tmScore);
            } catch (StructureException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        int index = this.getIndex(points, ThresholdConfig.TM_THRES);

        if (index > 0) {
            RepeatContent foundRepeat = new RepeatContent(end + 1, end + index + trsLen);
            repeat.getRepeats().add(foundRepeat);
            repeat.sortByPosition();
            extendRepeatRight(repeat, templateTRs);
        }


    }

    public int getIndex(List<Double> points, double tmTHRES) {
        int index = -1;
        double maxTM = -1;
        for (double a : points) {
            maxTM = Math.max(maxTM, a);
        }
        if (maxTM > tmTHRES) {
            for (int i = 0; i < points.size(); i++) {
                if (points.get(i) == maxTM) {
                    index = i;
                    break;
                }
            }

        } else
            index = -1;
        return index;
    }

    public void extendRepeatLeft(Repeat repeat, List<RepeatContent> templateTRs) {

        // doing something here
        int start = repeat.getStart();
        // int end = repeat.getEnd();
        // System.out.println(start);
        int trsLen = (int) repeat.getAvgLength();
        RepeatContent firstUnit = repeat.getRepeats().get(0);
        Atom[] atomTemp = Fragement.getFragementsofAtoms(atoms,
                firstUnit.getStart(), firstUnit.getEnd());

        List<Double> points = new ArrayList<Double>();
        for (int window = trsLen; window < start && window < trsLen * 2; window++) {
            try {
                Atom[] atomCompare = Fragement.getFragementsofAtoms(atoms,
                        start - window, start - 1);

                AFPChain afpChain = mutilAlign.pairAlign(atomTemp, atomCompare);
                double tmScore = afpChain.getTMScore();
                points.add(tmScore);
            } catch (StructureException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        int index = -1;
        if (trsLen > 25) {
            //first the first point
            index = this.getIndex(points, 0.5);

        } else {
            index = this.getIndex(points, 0.17);
        }
        if (index > 0) {
            RepeatContent foundRepeat = new RepeatContent(start - index - trsLen, start - 1);
            repeat.getRepeats().add(foundRepeat);
            repeat.sortByPosition();
            extendRepeatLeft(repeat, templateTRs);
        }


    }

    boolean checkPattern(AFPChain afpChain, int abPositionA1, int abPositionA2) {

        // boolean isSame = false;
        Map<Integer, Integer> map = mutilAlign.toAlignedPairs(afpChain,
                abPositionA1, abPositionA2);

        int match = 0;
        int index = 0;
        for (Entry<Integer, Integer> entry : map.entrySet()) {

            char charP1 = strSS.charAt(entry.getKey());
            char charP2 = strSS.charAt(entry.getValue());

            // if (charP1 != '-' && charP2 != '-') {
            // index++;
            //
            // if (charP1 == charP2)
            // match++;
            // }

            if (charP1 == charP2) {
                match++;
            }
            index++;

        }

        if ((double) match / index > 0.8)
            return true;
        else

            return false;

    }
}
