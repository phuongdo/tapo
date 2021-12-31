package org.cnrs.crbm.lib.repeats;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureTools;

/**
 * Fragment a protein
 *
 * @author phuongdo
 */
public class Fragement {
    public static Atom[] getFragementsofAtoms(Atom[] setAtoms, int posStart,
                                              int posEnd) {

        int size = posEnd - posStart + 1;
        Atom[] setAtoms1 = new Atom[size];
        // System.out.println(posStart + ":" + posEnd + " :" + size + ":atoms:"
        // + setAtoms.length);
        // clone
        for (int k = 0; k < size; k++) {
            // System.out.println(k);
            //setAtoms1[k] = (Atom)  setAtoms[posStart + k].clone();
            setAtoms1[k] = setAtoms[posStart + k];
        }
        return setAtoms1;
    }

    public static Atom[] getFragementsofAtomsWithClone(Atom[] setAtoms, int posStart,
                                                       int posEnd) {

        int size = posEnd - posStart + 1;
        Atom[] setAtoms1 = new Atom[size];
        // System.out.println(posStart + ":" + posEnd + " :" + size + ":atoms:"
        // + setAtoms.length);
        // clone
        for (int k = 0; k < size; k++) {
            // System.out.println(k);
            // setAtoms1[k] = (Atom) setAtoms[posStart + k].clone();
            setAtoms1[k] = setAtoms[posStart + k];
        }
        return StructureTools.cloneAtomArray(setAtoms1);

    }
//    public static Atom[] getAtomsWithSeqN(Atom[] atoms, int startNsq, int endNsq) {
//        List<Atom> subAtoms = new ArrayList<Atom>();
//        for (Atom atom : atoms) {
//            int seqNumber = atom.getGroup().getResidueNumber().getSeqNum();
//            if (seqNumber >= startNsq && seqNumber <= endNsq)
//                subAtoms.add(atom);
//        }
//        return subAtoms.toArray(new Atom[subAtoms.size()]);
//    }
}
