package org.cnrs.crbm.lib.align.unpublish;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.PDBFileReader;

public class AlignTest {

	public static Atom[] getAtomCBArray(Chain s) {
		String[] atomNames = { " CB " };
		return StructureTools.getAtomArray(s, atomNames);
	}

	public static void main(String[] args) {

		String pdb1 = "input/197_249.pdb";
		String pdb2 = "input/obj01.pdb";
		// String pdb1 = "input/1BYO.pdb";
		// String pdb2 = "input/1BYO.pdb";
		// String pdb1 = "input/1BPO.pdb";
		// String pdb2 = "input/4IRT.pdb";

		PDBFileReader pdbreader = new PDBFileReader();

		byte a = 1;

		Atom[] atomSet1 = null;
		Atom[] atomSet2 = null;
		try {
			Structure struc1 = pdbreader.getStructure(pdb1);
			Chain chain = struc1.getChainByPDB("A");
			atomSet1 = StructureTools.getAtomCAArray(chain);// ContactMapMain.getAtomCBArray(chain);
			Structure struc2 = pdbreader.getStructure(pdb2);
			// Chain chain2 = struc2.getChain(0);
			Chain chain2 = struc2.getChainByPDB("A");
			atomSet2 = StructureTools.getAtomCAArray(chain2);// ContactMapMain.getAtomCBArray(chain2);
		} catch (Exception e) {
			e.printStackTrace();
		}

		PairwiseStructureAlignment structureAlignment = new PairwiseStructureAlignment();
		AlignChain alignChain = structureAlignment.align(atomSet1, atomSet2);
		System.out.println(alignChain);

	}
}
