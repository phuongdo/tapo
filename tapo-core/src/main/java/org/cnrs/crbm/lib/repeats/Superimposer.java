package org.cnrs.crbm.lib.repeats;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;
import org.cnrs.crbm.lib.conf.ThresholdConfig;
import org.cnrs.crbm.lib.math.Threshold;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.module.ProVector;

import java.util.ArrayList;
import java.util.List;

public class Superimposer {
    /**
     * Just superimpose 2 array of Atom based on transformation and rotation of
     * their coordinates. This method is only affect on small size of two
     * structures. You should use CA atoms to superimpose.
     *
     * @param atomSet1
     * @param atomSet2
     * @throws StructureException
     */

    private static final double DIS_DIFF_THRESHOLD = 4.0;

    private boolean is_adjust = false;

    public double superimposeSimple(Atom[] atomSet1, Atom[] atomSet2)
            throws StructureException {

        if (atomSet1.length != atomSet2.length) {
            throw new StructureException(
                    "The two atom sets are not of same length!");
        }

        // clone before superimpose because of transformation and rotation

        Atom[] atomSet1Clone = new Atom[atomSet1.length];
        Atom[] atomSet2Clone = new Atom[atomSet2.length];
        for (int i = 0; i < atomSet1.length; i++) {
            atomSet1Clone[i] = (Atom) atomSet1[i].clone();
            atomSet2Clone[i] = (Atom) atomSet2[i].clone();
        }

        /**
         * START MODIFIDING
         */

        if (is_adjust) {
            double[] dists1 = this.getVectorDistances(atomSet1);
            double[] dists2 = this.getVectorDistances(atomSet2);
            int winsize = dists1.length;
            for (int i = 0; i < winsize; i++) {

                if (Math.abs(dists1[i] - dists2[i]) > DIS_DIFF_THRESHOLD) {

                    // Atom atom1_start = atomSet1Clone[i * 2];
                    // Atom atom1_end = atomSet1Clone[i * 2 + 1];
                    Atom atom2_start = atomSet2Clone[i * 2];
                    Atom atom2_end = atomSet2Clone[i * 2 + 1];
                    // get the vector unit
                    Point unit = this.getUnitVector(atom2_start, atom2_end);
                    // change the length
                    double x = atom2_start.getX() + dists1[i] * unit.getX();
                    double y = atom2_start.getY() + dists1[i] * unit.getY();
                    double z = atom2_start.getZ() + dists1[i] * unit.getZ();
                    atom2_end.setX(x);
                    atom2_end.setY(y);
                    atom2_end.setZ(z);
                    atomSet2Clone[i * 2 + 1] = atom2_end;
                }

            }
        }

        /**
         * END
         */

        SVDSuperimposer svds = new SVDSuperimposer(atomSet1Clone, atomSet2Clone);
        Matrix rotMatrix = svds.getRotation();
        Atom tranMatrix = svds.getTranslation();
        rotateShiftAtoms(atomSet2Clone, rotMatrix, tranMatrix);
        //return AFPChainScorerAlter.getTMScore(atomSet1Clone,atomSet2.clone(),atomSet1.length,atomSet2.length);
        return SVDSuperimposer.getRMS(atomSet1Clone, atomSet2Clone);

    }

    /**
     * a->b
     *
     * @param a
     * @param b
     * @return
     */
    public Point getUnitVector(Atom a, Atom b) {

        Point unit = new Point(a.getX() - b.getX(), a.getY() - b.getY(),
                a.getZ() - b.getZ());

        double norm = Math.sqrt(Math.pow(unit.getX(), 2)
                + Math.pow(unit.getY(), 2) + Math.pow(unit.getZ(), 2));

        unit.setX(unit.getX() / norm);
        unit.setY(unit.getY() / norm);
        unit.setZ(unit.getZ() / norm);

        return unit;
    }

    public void rotateShiftAtoms(Atom[] ca, Matrix rotMatrix, Atom tranMatrix) {

        for (int i = 0; i < ca.length; i++) {
            Atom c = (Atom) ca[i].clone();
            Calc.rotate(c, rotMatrix);
            Calc.shift(c, tranMatrix);
            ca[i] = c;
        }
    }

    private double[] getVectorDistances(Atom[] atoms) {
        int winsize = atoms.length / 2;
        double[] ref_dists = new double[winsize];
        for (int i = 0; i < winsize; i = i + 2) {
            Atom atom1 = atoms[i];
            Atom atom2 = atoms[i + 1];
            ref_dists[i] = Calc.getDistance(atom1, atom2);
        }
        return ref_dists;
    }

    public double superimposeAtoms(List<Atom> latom1, List<Atom> latom2)
            throws StructureException {
        Atom[] atomSet1 = latom1.toArray(new Atom[latom1.size()]);
        Atom[] atomSet2 = latom2.toArray(new Atom[latom2.size()]);

        return this.superimposeSimple(atomSet1, atomSet2);

    }


    /**
     * Hierarchical Protein Structure Superposition using both Secondary
     * Structure and Atomic Representations
     * <p/>
     * <pre>
     * Orientation Independent Scores:
     * S1 = S{|angle(i,k) - angle(p,Rlib)|}
     * S2 = S{|angle(i,j) - angle(p,q)|}
     * S3 = S{|angle(j,k) - angle(q, Rlib)|}
     * S4 = S{|distance(i,k) - distance(p,Rlib)|}
     * S5 = S{|length(k) - length(Rlib)|}
     * Orientation Dependent Scores:
     * S6 = S{angle(k,Rlib)}
     * S7 = S{distance(k,Rlib)}
     * </pre>
     *
     * @param
     * @throws Exception
     */

    public double compareListVector(List<ProVector> lvectors1,
                                    List<ProVector> lvectors2) throws Exception {

        double score = 0.0;
        if (lvectors1.size() != lvectors2.size())
            throw new Exception(" Two list of vectors must be the same size!!!");

        int size = lvectors1.size();

        try {


            for (int i = 0; i < size; i++) {

                double ith_score = 0.0;
                for (int j = 0; j < size; j++) {
                    // a pairs
                    if (j != i) {
                        ith_score += scoreTwoVector(lvectors1.get(i), lvectors1.get(j), lvectors2.get(i), lvectors2.get(j));

                    }
                }

                score += ith_score / (size - 1);


            }

//            for (int index = 0; index < size - 1; index++) {
//
//                ProVector i = lvectors1.get(index);
//                ProVector k = lvectors1.get(index + 1);
//
//                ProVector p = lvectors2.get(index);
//                ProVector r = lvectors2.get(index + 1);
//                score += scoreTwoVector(i, k, p, r);
//
//            }
        } catch (Exception ex) {
            // do something

            ex.printStackTrace();
        }
        return score / size;
    }


    private double scoreTwoVector(ProVector i, ProVector k, ProVector p,
                                  ProVector r) {

        double DIS_D0 = 2.0;
        double ANG_D0 = 10.0;

        Vector3D v1 = new Vector3D(i.getStart().getCoords());
        Vector3D v2 = new Vector3D(i.getEnd().getCoords());
        Vector3D v3 = new Vector3D(k.getStart().getCoords());
        Vector3D v4 = new Vector3D(k.getEnd().getCoords());
        Vector3D v12 = v2.subtract(v1);
        Vector3D v23 = v3.subtract(v2);
        Vector3D v34 = v4.subtract(v3);
        Vector3D totalV1 = v1.add(v2);
        Vector3D totalV2 = v3.add(v4);
        Vector3D mid_V1 = new Vector3D(totalV1.getX() / 2, totalV1.getY() / 2, totalV1.getZ() / 2);
        Vector3D mid_V2 = new Vector3D(totalV2.getX() / 2, totalV2.getY() / 2, totalV2.getZ() / 2);

        Vector3D v1_c = new Vector3D(p.getStart().getCoords());
        Vector3D v2_c = new Vector3D(p.getEnd().getCoords());
        Vector3D v3_c = new Vector3D(r.getStart().getCoords());
        Vector3D v4_c = new Vector3D(r.getEnd().getCoords());
        Vector3D v12_c = v2_c.subtract(v1_c);
        Vector3D v23_c = v3_c.subtract(v2_c);
        Vector3D v34_c = v4_c.subtract(v3_c);
        Vector3D totalV1_c = v1_c.add(v2_c);
        Vector3D totalV2_c = v3_c.add(v4_c);
        Vector3D mid_V1_c = new Vector3D(totalV1_c.getX() / 2, totalV1_c.getY() / 2, totalV1_c.getZ() / 2);
        Vector3D mid_V2_c = new Vector3D(totalV2_c.getX() / 2, totalV2_c.getY() / 2, totalV2_c.getZ() / 2);


//        if(v34_c.getNorm()==0){
//            System.out.println();
//        }

        double S1 = getScore_S(
                Math.abs(Math.toDegrees(Vector3D.angle(v12, v34))
                        - Math.toDegrees(Vector3D.angle(v12_c, v34_c))), 10, ANG_D0);
        double S2 = 1;
        if (v23.getNorm() != 0 && v23_c.getNorm() != 0) {
            S2 = getScore_S(
                    Math.abs(Math.toDegrees(Vector3D.angle(v12, v23))
                            - Math.toDegrees(Vector3D.angle(v12_c, v23_c))), 4, ANG_D0);
        }
        double S3 = 1;
        if (v23.getNorm() != 0 && v23_c.getNorm() != 0) {
            S3 = getScore_S(
                    Math.abs(Math.toDegrees(Vector3D.angle(v23, v34))
                            - Math.toDegrees(Vector3D.angle(v23_c, v34_c))), 4, ANG_D0);
        }


        /**
         * WRONG HERE, WE NEED TO RECALCULATE THIS DISTANCE
         * distance between two vectors is computed by averaging the distance between corresponding
         * start, end and middle point on the vector.
         */


        double distV12 = (Vector3D.distance(v1, v3) + Vector3D.distance(v2, v4) + Vector3D.distance(mid_V1, mid_V2)) / 3;
        double distV34 = (Vector3D.distance(v1_c, v3_c) + Vector3D.distance(v2_c, v4_c) + Vector3D.distance(mid_V1_c, mid_V2_c)) / 3;

        double S4 = getScore_S(
                Math.abs(distV12
                        - distV34), 2, DIS_D0);
        double S5 = getScore_S(Math.abs(v12.getNorm() - v12_c.getNorm()), 5, DIS_D0);
        double S6 = getScore_S(Math.abs(v34.getNorm() - v34_c.getNorm()), 5, DIS_D0);
        double S7 = 1;
        if (v23.getNorm() != 0 && v23_c.getNorm() != 0) {
            S7 = getScore_S(Math.abs(v23.getNorm() - v23_c.getNorm()), 5, DIS_D0);
        }

        return (S1 + S2 + S3 + S4 + S5 + S6 + S7) / 7;
        //return (S1 + S2 + S3 + S4) / 4;
    }

    private double getScore_S(double d, double M, double d0) {

        //return (2 * M / (1 + Math.pow(d / d0, 2)) - M);
        return (1 / (1 + Math.pow(d / d0, 2)));
    }

    /**
     * Avoid to use this method. Use : superimposeListVectors
     *
     * @param lvectors1
     * @param lvectors2
     * @return
     * @throws StructureException
     */
    @Deprecated
    public double superimposeVectors(List<ProVector> lvectors1,
                                     List<ProVector> lvectors2) throws StructureException {

        double rmsd = 10;
        if (lvectors1.size() == lvectors2.size()) {
            List<Atom> latom1 = this.convertToAtoms(lvectors1);
            List<Atom> latom2 = this.convertToAtoms(lvectors2);
            rmsd = this.superimposeAtoms(latom1, latom2);

        } else {

            int size1 = lvectors1.size();
            int size2 = lvectors2.size();
            if (size1 < size2) {
                List<Atom> latom1 = this.convertToAtoms(lvectors1);
                for (int i = 0; i < size2 - size1; i++) {
                    // find smallest rmsd
                    List<ProVector> l1 = new ArrayList<ProVector>();
                    for (int j = 0; j < size1; j++) {
                        l1.add(lvectors2.get(i + j));
                    }
                    double current_rmsd = this.superimposeAtoms(latom1,
                            this.convertToAtoms(l1));
                    if (current_rmsd < rmsd)
                        rmsd = current_rmsd;

                }

            } else {
                List<Atom> latom2 = this.convertToAtoms(lvectors2);
                for (int i = 0; i < size1 - size2; i++) {
                    // find smallest rmsd
                    List<ProVector> l2 = new ArrayList<ProVector>();
                    for (int j = 0; j < size2; j++) {
                        l2.add(lvectors1.get(i + j));
                    }
                    double current_rmsd = this.superimposeAtoms(latom2,
                            this.convertToAtoms(l2));
                    if (current_rmsd < rmsd)
                        rmsd = current_rmsd;

                }

            }

        }

        return rmsd;

    }

    private List<Atom> convertToAtoms(List<ProVector> lvectors) {
        List<Atom> latom1 = new ArrayList<Atom>();
        for (ProVector v : lvectors) {
            latom1.add(v.getStart());
            latom1.add(v.getEnd());
        }
        return latom1;
    }

    public SuperimposerOutput compareTwoStructures(Atom[] atomSet1, Atom[] atomSet2)
            throws StructureException {

        SuperimposerOutput output = new SuperimposerOutput();
        try {
            MutilAlign mutilAlign = new MutilAlign();
            int noRes1 = atomSet1.length;
            int noRes2 = atomSet2.length;
            if (noRes1 > 10) {
                AFPChain afpChain = null;
                afpChain = mutilAlign.pairAlign(atomSet1, atomSet2);
                double alnLeng = (double) afpChain.getOptLength();
                double localMaxRMSD = Threshold.functRMSD(afpChain
                        .getOptLength());
                double aL = Threshold.functaL(Math.max(noRes1, noRes2));
                double aS = Threshold.functaS(Math.min(noRes1, noRes2));
                double gapsM = Threshold.functGapsRes(alnLeng);
                double tmScoreMin = Threshold.funcTmScore(alnLeng);
                int gaps = afpChain.getGapLen();
                double tmScore = afpChain.getTMScore();
                //System.out.println(tmScore);
                // tmScore = 1;
                int Lmin = Math.min(noRes1, noRes2);

                //System.out.println(tmScore);
                output.setTmScore(tmScore);

                /**
                 * TESTING ???????
                 */
                if (Lmin < 60) {
                    if (tmScore >= ThresholdConfig.TM_THRES && alnLeng > 0.7 * Lmin) {
                        output.setIsTRs(true);
                    } else
                        output.setIsTRs(false);
                } else {
                    /**
                     * Default TM for very long repeat like domain.
                     */
                    if (tmScore >= 0.5) output.setIsTRs(true);
                    else output.setIsTRs(false);
                }


            }
        } catch (Exception e) {
            // System.out.println("set 1: " + atomSet1.length + "  set 2: "
            // + atomSet2.length);
            e.printStackTrace();
        }
        return output;

    }


    /**
     * Improve the accuracy of superimpose algorithm by considering both length
     * of atoms and rmsd threshold.
     *
     * @param atomSet1 reference atom set
     * @param atomSet2 comparing atom set
     * @return
     * @throws StructureException
     */
    public boolean isTheSameStructure(Atom[] atomSet1, Atom[] atomSet2)
            throws StructureException {
        try {
            MutilAlign mutilAlign = new MutilAlign();
            int noRes1 = atomSet1.length;
            int noRes2 = atomSet2.length;
            if (noRes1 > 10) {
                AFPChain afpChain = null;


                afpChain = mutilAlign.pairAlign(atomSet1, atomSet2);

                double alnLeng = (double) afpChain.getOptLength();
                double localMaxRMSD = Threshold.functRMSD(afpChain
                        .getOptLength());
                double aL = Threshold.functaL(Math.max(noRes1, noRes2));
                double aS = Threshold.functaS(Math.min(noRes1, noRes2));
                double gapsM = Threshold.functGapsRes(alnLeng);
                double tmScoreMin = Threshold.funcTmScore(alnLeng);
                int gaps = afpChain.getGapLen();
                double tmScore = afpChain.getTMScore();
                //System.out.println(tmScore);
                // tmScore = 1;
                int Lmin = Math.min(noRes1, noRes2);

                //System.out.println(tmScore);

                // dommain
                if (Lmin < 60) {
                    if (tmScore >= ThresholdConfig.TM_THRES) return true;
                    else return false;
                } else {
                    /**
                     * Default TM for very long repeat like domain.
                     */
                    if (tmScore >= 0.5) return true;
                    else return false;
                }


//                int minL = Math.min(noRes1, noRes2);
//                int maxL = Math.max(noRes1, noRes2);
//                if (minL > 30) {
//                    if (tmScore >= 0.5)
//                        return true;
//                    else return false;
//                } else {
//
//                    // perfect matched
//
//                    if (tmScore > 0.17 && afpChain.getChainRmsd() <= localMaxRMSD && gaps <= gapsM
//                            && alnLeng >= aS && alnLeng >= aL) {
//                        return true;
//
//                    } else return false;
//
//
////                    if (alnLeng == minL && maxL == minL)
////                        return true;
////                    else return false;
//
//
//                }


//                if (afpChain.getChainRmsd() <= localMaxRMSD && gaps <= gapsM
//                        && alnLeng >= aS && alnLeng >= aL && tmScore > tmScoreMin) {
//
////					if (Math.min(noRes1, noRes2) < 30)
////						return true;
////					else {
////						if (tmScore > 0.5)
////							return true;
////						else
////							return false;
////					}
//                    return true;
//                } else
//                    return false;

                // end < localMaxRMD
            } else
                return false;
        } catch (Exception e) {
            // System.out.println("set 1: " + atomSet1.length + "  set 2: "
            // + atomSet2.length);
            e.printStackTrace();
        }
        return false;

    }

    public boolean isTheSameContactMap(Atom[] atomSet1, int startSet1,
                                       Atom[] atomSet2, int startSet2, double[] cmHistos)
            throws StructureException {
        try {
            MutilAlign mutilAlign = new MutilAlign();
            AFPChain afpChain = mutilAlign.pairAlign(atomSet1, atomSet2);
            double rmsdOfContactMap = mutilAlign.alignContactMap(afpChain,
                    startSet1, startSet2, cmHistos);
            double localMaxRMSDContact = 0.6;// fix later
            if (rmsdOfContactMap < localMaxRMSDContact) {
                return true;
            }
            return false;
        } catch (Exception e) {
            // System.out.println("set 1: " + atomSet1.length + "  set 2: "
            // + atomSet2.length);
            e.printStackTrace();
        }
        return false;

    }
}
