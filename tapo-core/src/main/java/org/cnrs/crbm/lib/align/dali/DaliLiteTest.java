package org.cnrs.crbm.lib.align.dali;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.util.Map;



public class DaliLiteTest {

	public static void main(String[] args) throws Exception {

		DaliLiteTest test = new DaliLiteTest();
		// test.test();
		// System.exit(0);
		String pdb1 = args[0];
		String pdb2 = args[1];

		JavaWapperDali daliLite = new JavaWapperDali();

		PDBFileReader pdbreader = new PDBFileReader();
		Atom[] atomSet1 = null;
		Atom[] atomSet2 = null;
		try {

			Chain struc1 = pdbreader.getStructure(pdb1).getChain(0);
			Chain struc2 = pdbreader.getStructure(pdb2).getChain(0);
			atomSet1 = StructureTools.getAtomCAArray(struc1);
			atomSet2 = StructureTools.getAtomCAArray(struc2);
		} catch (Exception e) {
			e.printStackTrace();
		}

		AlignChainDali alignChain = daliLite.align(atomSet1, atomSet2);

		System.out.println(alignChain.getChainLength());
		System.out.println(alignChain.getRmsd());
		System.out.println(alignChain.getGaps());
		Map<Integer, Integer> aligns = alignChain.getAlignPairs();
		for (Map.Entry<Integer, Integer> entry : aligns.entrySet()) {
			System.out.println("p1 = " + entry.getKey() + ", p2 = "
					+ entry.getValue());
		}

	}

	public void test() {

		DaliIO daliIO = new DaliIO();
		AlignChainDali alignChain = daliIO.readOutput("input/summary.txt");
		System.out.println(alignChain.getChainLength());
		System.out.println(alignChain.getRmsd());
		System.out.println(alignChain.getGaps());
		Map<Integer, Integer> aligns = alignChain.getAlignPairs();
		for (Map.Entry<Integer, Integer> entry : aligns.entrySet()) {
			System.out.println("p1 = " + entry.getKey() + ", p2 = "
					+ entry.getValue());
		}
	}
}
