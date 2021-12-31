package org.cnrs.crbm.lib.align.dali;

import java.util.Map;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.utils.ExecUtils;

public class JavaWapperDali {
	DaliIO daliIO = new DaliIO();

	private AlignChainDali callDaliLite(String pdbFile1, String pdbFile2) {

		pdbFile1 = Dir.DALI_TMP_DIR + "/" + pdbFile1 + ".pdb";
		pdbFile2 = Dir.DALI_TMP_DIR + "/" + pdbFile2 + ".pdb";

		// remove lock first
		ExecUtils.execShellCmdLinux("rm -rf " + Dir.HOME_DIR + "/dali.lock");
		ExecUtils.execShellCmdLinux("rm -rf " + Dir.HOME_DIR + "/summary.txt");
		// ExecUtils.execShellCmd2("rm -rf " + Dir.DALI_TMP_DIR + "/dali.lock");
		// ExecUtils.execShellCmd2("rm -rf " + Dir.DALI_TMP_DIR +
		// "/summary.txt");

		String shellCommand = Dir.DALI_DIR + "/DaliLite -pairwise " + pdbFile1
				+ " " + pdbFile2;

		// run dali
		String strCmdOut = ExecUtils.execShellCmdLinux(shellCommand);

		// read output
		String outputFile = Dir.HOME_DIR + "/summary.txt";
		AlignChainDali alignChain = daliIO.readOutput(outputFile);
		return alignChain;

	}

	public AlignChainDali align(Atom[] atomSet1, Atom[] atomSet2)
			throws StructureException {

		// save atom to pdbFile
		Map<Integer, Integer> mapStruct1 = daliIO.writePDB(atomSet1, "tmp1");
		Map<Integer, Integer> mapStruct2 = daliIO.writePDB(atomSet2, "tmp2");

		AlignChainDali alignChain = this.callDaliLite("tmp1", "tmp2");

		// System.out.println("----Dali call----");
		// System.out.println(alignChain.getChainLength());
		// System.out.println(alignChain.getRmsd());
		// System.out.println(alignChain.getGaps());
		// Map<Integer, Integer> aligns = alignChain.getAlignPairs();
		// for (Map.Entry<Integer, Integer> entry : aligns.entrySet()) {
		// System.out.println("p1 = " + entry.getKey() + ", p2 = "
		// + entry.getValue());
		// }
		// System.out.println("----end----");
		return alignChain;

	}
}
