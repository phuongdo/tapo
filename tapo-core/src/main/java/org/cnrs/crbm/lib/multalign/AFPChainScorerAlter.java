package org.cnrs.crbm.lib.multalign;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Created by pdoviet on 11/5/2014.
 */
public class AFPChainScorerAlter {


    public static double getTMScore(AFPChain align, Atom[] ca1, Atom[] ca2) throws StructureException {
        if (align.getNrEQR() == 0)
            return -1;


        // Create new arrays for the subset of atoms in the alignment.
        Atom[] ca1aligned = new Atom[align.getOptLength()];
        Atom[] ca2aligned = new Atom[align.getOptLength()];
        int pos = 0;
        int[] blockLens = align.getOptLen();
        int[][][] optAln = align.getOptAln();
        assert (align.getBlockNum() <= optAln.length);

        for (int block = 0; block < align.getBlockNum(); block++) {

            if (!(blockLens[block] <= optAln[block][0].length)) {
                System.err.println("AFPChainScorer getTMScore: errors reconstructing alignment block [" + block + "]. Length is " + blockLens[block] + " but should be <=" + optAln[block][0].length);
            }

            for (int i = 0; i < blockLens[block]; i++) {
                int pos1 = optAln[block][0][i];
                int pos2 = optAln[block][1][i];
//				if ( pos1 >= ca1.length) {
//					System.err.println("getTMScore: pos1 " + pos1 + " > ca1.length " + ca1.length + " " + pairAlign.getName1() + ":" +pairAlign.getName2());
//					System.err.println(ca1[0].getGroup().getChain().getChainID());
//					System.err.println(ca1[0].getGroup().getChain().getParent().getPDBCode());
//
//				}
//				if ( pos2 >= ca2.length) {
//					System.err.println("getTMScore: pos2 " + pos2 + " > ca2.length " + ca2.length + " " + pairAlign.getName2() + ":" + pairAlign.getName1());
//					System.err.println(ca2[0].getGroup().getChain().getChainID());
//					System.err.println(ca2[0].getGroup().getChain().getParent().getPDBCode());
//				}
                Atom a1 = ca1[pos1];
                Atom a2 = (Atom) ca2[pos2].clone();

                ca1aligned[pos] = a1;
                ca2aligned[pos] = a2;
                pos++;
            }
        }

        // this can happen when we load an old XML serialization which did not support modern ChemComp representation of modified residues.
        if (pos != align.getOptLength()) {
            System.err.println("AFPChainScorer getTMScore: Problems reconstructing alignment! nr of loaded atoms is " + pos + " but should be " + align.getOptLength());
            // we need to resize the array, because we allocated too many atoms earlier on.
            ca1aligned = (Atom[]) resizeArray(ca1aligned, pos);
            ca2aligned = (Atom[]) resizeArray(ca2aligned, pos);
        }
        //Superimpose
        SVDSuperimposer svd = new SVDSuperimposer(ca1aligned, ca2aligned);
        Matrix matrix = svd.getRotation();
        Atom shift = svd.getTranslation();

        for (Atom a : ca2aligned) {
            Calc.rotate(a, matrix);
            Calc.shift(a, shift);
        }

        return getTMScore(ca1aligned, ca2aligned, ca1.length, ca2.length);

    }

    /**
     * Calculate the TM-Score for the superposition.
     * <p/>
     * <em>Normalizes by the <strong>minimum</strong>-length structure (that is, {@code min\{len1,len2\}}).</em>
     * <p/>
     * Atom sets must be pre-rotated.
     * <p/>
     * <p>Citation:<br/>
     * <i>Zhang Y and Skolnick J (2004). "Scoring function for automated assessment
     * of protein structure template quality". Proteins 57: 702 - 710.</i>
     *
     * @param atomSet1 atom array 1
     * @param atomSet2 atom array 2
     * @param len1     The full length of the protein supplying atomSet1
     * @param len2     The full length of the protein supplying atomSet2
     * @return The TM-Score
     * @throws StructureException
     */
    public static double getTMScore(Atom[] atomSet1, Atom[] atomSet2, int len1, int len2) throws StructureException {
        if (atomSet1.length != atomSet2.length) {
            throw new StructureException("The two atom sets are not of same length!");
        }
        if (atomSet1.length > len1) {
            throw new StructureException("len1 must be greater or equal to the alignment length!");
        }
        if (atomSet2.length > len2) {
            throw new StructureException("len2 must be greater or equal to the alignment length!");
        }

        int Lmin = Math.min(len1, len2);
        int Laln = atomSet1.length;

        double D0_MIN = 0.5;
        double d0 = 1.24 * Math.cbrt(Lmin - 15.) - 1.8;
        //if (d0 < D0_MIN) d0 = D0_MIN;
        d0 = 1.24 * Math.cbrt(Lmin - 10.) - 1.8;
//        if (Lmin < 25) {
//            d0 = 1.24 * Math.cbrt(Lmin - 10.) - 1.8;
//            //d0 = 0.5;
//
//        }
        double d0sq = d0 * d0;

        double sum = 0;
        for (int i = 0; i < Laln; i++) {
            double d = Calc.getDistance(atomSet1[i], atomSet2[i]);
            sum += 1. / (1 + d * d / d0sq);
        }

        return sum / Lmin;
    }

    /**
     * Reallocates an array with a new size, and copies the contents
     * of the old array to the new array.
     *
     * @param oldArray the old array, to be reallocated.
     * @param newSize  the new array size.
     * @return A new array with the same contents.
     */
    private static Object resizeArray(Object oldArray, int newSize) {
        int oldSize = java.lang.reflect.Array.getLength(oldArray);
        @SuppressWarnings("rawtypes")
        Class elementType = oldArray.getClass().getComponentType();
        Object newArray = java.lang.reflect.Array.newInstance(
                elementType, newSize);
        int preserveLength = Math.min(oldSize, newSize);
        if (preserveLength > 0)
            System.arraycopy(oldArray, 0, newArray, 0, preserveLength);
        return newArray;
    }


}
