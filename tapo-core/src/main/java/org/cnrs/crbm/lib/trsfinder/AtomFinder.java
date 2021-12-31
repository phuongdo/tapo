package org.cnrs.crbm.lib.trsfinder;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Deprecated
public class AtomFinder {
    static Logger logger = LoggerFactory.getLogger(AtomFinder.class);

    public static List<Repeat> getRepeat(Features feature, int repeatLeng) {
        return getPureRMSDFinder(feature, repeatLeng, 0,
                feature.getAtoms().length - 1);

    }

    public static List<Repeat> getRepeatRange(Features features,
                                              int repeatLeng, int start, int end) {
        return getPureRMSDFinder(features, repeatLeng, start, end);

    }

    public static List<Repeat> getPureRMSDFinder(Features features,
                                                 int winsize, int startFragment, int endFragment) {

        Superimposer superimposer = new Superimposer();
        List<Repeat> pures = new ArrayList<Repeat>();
        Atom[] atoms = features.getAtoms();
        String strSS = features.getStrSS();
        // pures.addAll(AtomFinder.getRepeat(features, winsize));
        int lastTRs = -1;
        Atom[] as1 = null;
        Atom[] as2 = null;
        for (int i = startFragment; i < endFragment - winsize; i++) {

            if (lastTRs >= 0 && lastTRs + winsize > i) {
                i = lastTRs + winsize - 1;
                continue;

            }
            org.cnrs.crbm.lib.trsfinder.Repeat aRepeat = new org.cnrs.crbm.lib.trsfinder.Repeat();
            as1 = Fragement.getFragementsofAtoms(atoms, i, i + winsize - 1);
            aRepeat.getRepeats().add(new RepeatContent(i, i + winsize - 1));
            String pattern1 = VectorShape.getSSPattern(strSS
                    .substring(i, i + winsize - 1));
            int lastFound = i + winsize;
            for (int j = i + winsize; j < endFragment - winsize; ) {

                as2 = Fragement.getFragementsofAtoms(atoms, j, j + winsize - 1);
                // compare
                try {
                    String pattern2 = VectorShape.getSSPattern(strSS
                            .substring(j, j + winsize - 1));
                    // check pattern
                    boolean check = true;
                    if (as1.length < 25) {
                        check = pattern1.equals(pattern2) && (pattern1.length() > 1 || pattern1.length() == 0);
                    } else if (as1.length < 40) {
                        if (pattern1.length() > 1 && pattern1.length() <= 3)
                            check = pattern1.equals(pattern2);
                        else if (pattern1.length() == 1 && pattern1.equals(pattern2))
                            check = false;
                        else
                            check = true;
                    } else if (pattern1.length() == 1 && pattern1.equals(pattern2)) {
                        check = false;
                    }

                    if (check && superimposer.isTheSameStructure(as1, as2)) {
                        aRepeat.getRepeats().add(
                                new RepeatContent(j, j + winsize - 1));

                        j = j + winsize - 1;
                        lastFound = j;

                    } else {
                        j++;
                        if (j - lastFound > (double) 0.2 * winsize)// if it's
                            // too far
                            // from last
                            // point, break
                            break;
                    }
                } catch (StructureException e) {

                    logger.error(e.getMessage() + " :" + features.getPdbCode());
                    // e.printStackTrace();
                    j++;
                }

            }

            if (aRepeat.getRepeats().size() >= 2) {
                pures.add(aRepeat);
                lastTRs = i;
            }

        }

        return pures;

    }

    public static List<Repeat> getRepeatOld(Features feature, int repeatLeng) {
        Superimposer superimpose = new Superimposer();
        Atom[] atoms = feature.getAtoms();
        List<ProVector> vectors = ProVector.toSecondaryVector(feature
                .getVectors());
        // String strSS = feature.getStrSS();
        List<Repeat> repeats = new ArrayList<Repeat>();
        for (int i = 0; i < vectors.size() - 1; i++) {
            ProVector v1 = vectors.get(i);
            ProVector v2 = vectors.get(i + 1);
            // get the middle of secondary structure element
            int seedPos = (v1.getPosEnd() + v2.getPosStart()) / 2;
            int nextPos = -1;
            String seedPattern = "";
            // find nextPos
            int k = -1;
            for (int j = i + 1; j < vectors.size() - 1; j++) {
                seedPattern += vectors.get(j).getType();
                nextPos = (vectors.get(j).getPosEnd() + vectors.get(j + 1)
                        .getPosStart()) / 2;
                if (Math.abs(nextPos - seedPos + 1 - repeatLeng) < 5) {
                    // found nextPos
                    k = j;
                    break;
                }
                if (nextPos - seedPos + 1 - repeatLeng > 5)
                    break;
            }

            if (nextPos == -1 || k == -1 || seedPattern.length() <= 1)
                continue;
            try {
                // supper impose with pattern??
                Atom[] seed = Fragement.getFragementsofAtoms(atoms, seedPos,
                        nextPos - 1);
                int nrSecondaryElement = seedPattern.length();

                if (k + nrSecondaryElement + 1 > vectors.size() - 1)
                    continue;
                int nextOfNextPos = (vectors.get(k + nrSecondaryElement)
                        .getPosEnd() + vectors.get(k + nrSecondaryElement + 1)
                        .getPosStart()) / 2;

                String comparePattern = "";
                for (int n = 0; n < nrSecondaryElement; n++) {
                    comparePattern += vectors.get(k + 1 + n).getType();
                }

                // System.out.println(nextPos + "-" + nextOfNextPos);
                Atom[] as2 = Fragement.getFragementsofAtoms(atoms, nextPos,
                        nextOfNextPos);

                if (seedPattern.equals(comparePattern)
                        && superimpose.isTheSameStructure(seed, as2)) {
                    Repeat a_repeat = new Repeat();
                    a_repeat.getRepeats().add(
                            new RepeatContent(seedPos, nextPos - 1));
                    a_repeat.getRepeats().add(
                            new RepeatContent(nextPos, nextOfNextPos));

                    // 2014/08/13
                    // after finding TRs for seed, we try to find more repeats
                    // based on these tandem repeats

                    for (int t = nextOfNextPos + 1; t < atoms.length
                            - repeatLeng; t = t + repeatLeng) {

                        Atom[] comparingAtoms = Fragement.getFragementsofAtoms(
                                atoms, t, t + repeatLeng - 1);
                        if (superimpose
                                .isTheSameStructure(seed, comparingAtoms)) {

                            a_repeat.getRepeats().add(
                                    new RepeatContent(t, t + repeatLeng - 1));

                        } else
                            // if it's two far, break
                            break;

                    }

                    repeats.add(a_repeat);

                }
            } catch (StructureException e) {
                // TODO Auto-generated catch block
                // e.printStackTrace();
            }
        }
        return repeats;
    }

}
