package org.cnrs.crbm.lib.sadb;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map.Entry;


import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.io.GzipIO;

import org.cnrs.crbm.lib.utils.PdbTools;

/**
 * The goal of defining a structural alphabet is to code a 3D structure fragment
 * of protein backbones and is to represent a 3D protein structure by a serial
 * of structural alphabets. An alphabet represents pattern profiles of the
 * backbone fragments (five residues long) derived from the pair database,
 * therefore, a protein structure of L residues is described by a structural
 * alphabet sequence of L-4 alphabets. We developed a nearest-neighbor
 * clustering (NNC) algorithm to cluster 225523 3D-protein fragments into 23
 * groups, which are represented by respective structural alphabets. We found
 * that these 23 structural alphabets can represent the profiles of most of the
 * 3D fragments and be roughly divided into five categories: Helix alphabet (A,
 * Y, B, C, and D), helix-like alphabet (G, I, and L), strand alphabet (E, F,
 * and H), strand-like alphabet (K and N), and others. The 3D sharps of
 * representation segments in the same category are similar. For example, the
 * sharps of 3D segments in the helix alphabets are similar and the ones of
 * strand alphabets are also similar. These 3D-fragment sharps and structural
 * alphabets are shown in the following Figure.
 * 
 * http://3d-blast.life.nctu.edu.tw/sa.php J.-M. Yang and C.-H. Tung
 * "Protein structure database search and evolutionary trclassification," Nucleic
 * Acids Research, vol. 34, pp. 3646-3659, 2006.
 * 
 * @author phuongdo
 * 
 */
public class StructAlphabet {

	public StructAlphabet() {
	}

	String proCode;
	String pdbSAFile;

	public StructAlphabet(String proCode) {

		this.pdbSAFile = Dir.SADB_LOCAL + "/"
				+ proCode.substring(1, 3).toLowerCase() + "/" + proCode + ".sa";

		try {

			File f = new File(this.pdbSAFile);
			if (!f.exists()) {

				insertDB(proCode);

			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void insertDB(String pdb) {

		String fileDir = Dir.SADB_LOCAL + "/"
				+ pdb.substring(1, 3).toLowerCase();
		File theDir = new File(fileDir);
		if (!theDir.exists()) {
			// make a directory before save file
			theDir.mkdir();
		}

		// extract file

		String pdbFileGz = Dir.PDB_LOCAL + "/" + pdb.substring(1, 3) + "/pdb"
				+ pdb + ".ent.gz";
		String tmpOutFilePDB = Dir.TMP_DIR + "/" + pdb + ".pdb";

		GzipIO.extract(pdbFileGz, tmpOutFilePDB);

		// save file here after running DSSP programme
		String outputFile = fileDir + "/" + pdb + ".sa";

		Structure structure = PdbTools.getStructureFromLocalPdb(pdb);

		for (Chain chain : structure.getChains()) {

			// this.calBashshell(tmpOutFilePDB, chain.getChainID(), outputFile);
		}
		// delete tmp file

		new File(tmpOutFilePDB).delete();

	}

	// private void generateSAfile() {
	//
	// Structure structure = PdbTools.getStructureFromLocalPdb(proCode);
	//
	// for (Chain chain : structure.getChains()) {
	//
	// // System.out.println(chain.getChainID());
	// // generate cmd
	// // this.calBashshell(pdbFile, chain.getChainID());
	//
	// }
	//
	// }

	// private void calBashshell(String pdbFile, String pdbChain, String
	// outFile) {
	//
	// String shellCommand = Dir.SADB_EXECUTABLE + " -sq_append " + pdbFile
	// + " " + pdbChain + " -o " + outFile;
	// System.out.println(shellCommand);
	// String strCmdOut = ExecUtils.execShellCmdLinux(shellCommand);
	//
	// }

	/**
	 * 
	 * @return the alphabet string of one pdb structure
	 * @throws Exception
	 */
	public String getAlphabets(String pdbCode, String chainCode)
			throws Exception {
		// Try with the FastaReaderHelper
		String strAlphabets = "";
		LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
				.readFastaProteinSequence(new File(this.pdbSAFile));
		for (Entry<String, ProteinSequence> entry : a.entrySet()) {
			if (entry.getValue().getOriginalHeader()
					.contains(pdbCode + ".pdb_" + chainCode)) {
				strAlphabets = entry.getValue().getSequenceAsString();
				break;
			}
		}
		if (strAlphabets.length() == 0)
			throw new Exception("alphabets structure of " + pdbCode
					+ " is not found!!!");
		return this.fixAlphabetProblems(strAlphabets);
	}

	final static char SA_coding_rule[] = "EEFFFFZZZZZZZZZZZZZZZZZZZZZZZZZZNNHHEEEFFFKKZZZZZZZZZZZZZZZZZZZZZXXXXNHHHHHKFFKKKKKZZZZZZZZZZZZZZZZZXXXXHHHHHNNNKKKKKKKKKZZZZZZZZZZZZXZXXXXXHHHHNPTTTNNKKKKKKKXKXXXXZXXXXXXXXXXXXNNNPPPPTTTTNNNTXXXXXXXXXXXXXXXXXXXXXRRPPPPPPTTTTTTTTTVVVVVXMXMMMMXXXRRRRRPPPPPPPPTTTTTTVVVVVVVMMMMMMQQQQQQRRQPPPPSSSSTTTTTTTVVVVVVVVMGMGGQQQQQQQQRQQSSSSSSSSSSVVVVVVVVVVMGGGGQQQQQQQQZZZZZSSSSSSSSSSWWWVVVVMDDBDLQQQQQRZZZZZZZZZSSSSSSSSWWWWWWWLDACIILLQRRZZZZZZZZZZZSSSSSSWWWWWWLLLLLIILLLRZZZZZZZZZZZZZZZSSSSSWWLWLLLLLIILZZZZZZZZZZZZZZZZZZZZZZZZZZZZZLZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
			.toCharArray();

	public String getSA(StrucState[] strucStates) {
		char line[] = new char[255];
		char line_chain = '-';
		char ftmp[] = new char[7];
		char[][] coding = new char[18][36];
		StringBuffer sabuffer = new StringBuffer();
		// String SA_string = "";
		char[] output = new char[1];
		double kappa = 0, kappa_tmp = 0, alpha = 0, alpha_tmp = 0;

		int k, i, j;
		k = 0;
		for (i = 0; i < 18; i++) {
			for (j = 0; j < 36; j++) {
				coding[i][j] = SA_coding_rule[k];
				k++;
			}
		}

		for (StrucState ss : strucStates) {

			kappa = ss.getKappa();
			alpha = ss.getAlpha();
			if (kappa == 180.0)
				kappa = 179.0;
			if (alpha == 180.0)
				alpha = 179.0;
			kappa = Math.floor(kappa / 10);
			alpha = Math.floor((alpha + 180) / 10);

			if (line[13] == '!')
				continue;

			if (kappa > 18 || alpha > 36) {
				// continue;
				output[0] = 'X';
			} else {
				output[0] = coding[(int) kappa][(int) alpha];
			}
			if (output[0] == 'A')
				if (kappa_tmp <= 114 && alpha_tmp > 46)
					output[0] = 'Y';

			sabuffer.append(output[0]);

		}

		return sabuffer.toString();
	}

	/**
	 * Because of generating alphabets structure, both 2 amino axit of start and
	 * end of a protein do not have their structure.
	 * 
	 * @param strAlphabets
	 * @return
	 */
	private String fixAlphabetProblems(String strAlphabets) {
		String nStr = strAlphabets.substring(2, strAlphabets.length() - 2);
		nStr = nStr.charAt(0) + "" + nStr.charAt(0) + "" + nStr + ""
				+ nStr.charAt(nStr.length() - 1) + ""
				+ nStr.charAt(nStr.length() - 1);
		return nStr;
	}

	public static void main(String[] args) throws Exception {
		String pdbId = "3lyc";
		String chainId = "A";
		if (args.length > 0) {
			pdbId = args[0];
			chainId = args[1];
		}

		Structure structure = PdbTools.getStructureFromLocalPdb(pdbId);
		Atom[] atoms = PdbTools
				.getAtomCAArray(structure.getChainByPDB(chainId));
		DSSP dssp = new DSSP(pdbId);
		StrucState[] strucStates = dssp.getStrucState(chainId, atoms);
		StructAlphabet structAlphabet = new StructAlphabet();
		System.out.println(structAlphabet.getSA(strucStates));
	}
}
