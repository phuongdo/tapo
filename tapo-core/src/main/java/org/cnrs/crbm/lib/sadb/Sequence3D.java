package org.cnrs.crbm.lib.sadb;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.utils.PdbTools;


/**
 * <pre>
 * Conformational Sequence and their corresponding torsion angles
 * 
 *  e (or k if inside)
 *     	phi +AD0- 40 to 160 and psi +AD0- -180:-40
 *         phi +AD0- 40 to 160 and psi +AD0- 145:180
 * 	
 *  l (or i if inside)
 * 	phi +AD0- 30 to 100 and psi +AD0- -20 to 90
 *       
 *  b (or v if inside)
 *         phi +AD0- -180 to -90 and psi +AD0-   90 to 180
 *         phi +AD0- -180 to -90 and psi +AD0- -180 to-165
 * 	phi +AD0-  160 to 180 and psi +AD0-   90 to 180
 *         phi +AD0-  160 to 180 and psi +AD0- -180 to-165
 *        
 *  p (or s if inside)
 * 	phi +AD0- -90 to -30 and psi +AD0-   90 to 180
 *         phi +AD0- -90 to -30 and psi +AD0- -180 to-165
 * 	
 *  d (or c if inside)
 * 	
 *         phi +AD0- -180 to -50 and psi +AD0-   20 to 90
 *         this is c if inside
 * 
 *  g (or h if inside)
 * 	phi +AD0- -160 to -50 and psi +AD0-   -25 to 20
 *                 
 *  a (or f if inside)
 *  	phi +AD0- -140 to -20 and psi +AD0- -80 to -25
 * 	this is f if inside
 * 
 *  x  - any other conformation
 *  
 *  yy - for two residues with peptide group having omega not equal to 180+-30,
 * </pre>
 * 
 * @author phuongdo
 * 
 */
public class Sequence3D {

	public static void main(String[] args) throws Exception {
		String pdbCode = "4mtl";
		String pdbChain = "A";

		if (args.length > 0) {
			pdbCode = args[0];
			pdbChain = args[1];
			// winsize = Integer.parseInt(args[2]);
		}

		// String pdbFile = PdbTools.downloadPDB(pdbCode);
		// Structure structure = PdbTools.getStructureFromFile(pdbFile);
		//

		Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
		Atom[] atoms = PdbTools.getAtomCAArray(structure
				.getChainByPDB(pdbChain));

		DSSP dssp = new DSSP(pdbCode);
		Sequence3D sequence3D = new Sequence3D();
		String strAccs = dssp.getRelativeACC(pdbChain, atoms);
		StrucState[] strucStates = dssp.getStrucState(pdbChain, atoms);
		System.out.println(atoms.length);
		System.out.println(strAccs.length());

		String sa = sequence3D.getSA(strucStates, strAccs);
		System.out.println(sa);
		System.out.println("converting no inside or outside");
		System.out.println(sequence3D.getSAV2NoInOutSide(sa));
		System.out.println("converting no inside or outside");
		System.out.println(sequence3D.getSAV3(sa));

	}

	/**
	 * Generate structural alphabets record
	 * 
	 * @param pdbCode
	 * @param pdbChain
	 * @return
	 */
	public String getSARecord(String pdbCode, String pdbChain) {

		StringBuffer buffer = new StringBuffer();
		try {
			Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
			Atom[] atoms;

			atoms = PdbTools.getAtomCAArray(structure.getChainByPDB(pdbChain));

			DSSP dssp = new DSSP(pdbCode);
			Sequence3D sequence3D = new Sequence3D();
			String strAccs = dssp.getRelativeACC(pdbChain, atoms);
			StrucState[] strucStates = dssp.getStrucState(pdbChain, atoms);

			buffer.append(">" + pdbCode + "_" + pdbChain
					+ "|PDBID|CHAIN|ALPHABETS SEQUENCE\n");
			buffer.append(sequence3D.getSA(strucStates, strAccs) + "\n");
		} catch (Exception e) {

			// System.out.println(pdbCode + "_" + pdbChain);
			e.printStackTrace();
			// clear
			buffer.delete(0, buffer.length());
		}
		return buffer.toString();

	}

	/**
	 * e l d p b a bb aa
	 * 
	 * @param sa
	 * @return
	 */
	public String getSAV2NoInOutSide(String sa) {

		sa = sa.replaceAll("k", "e");
		sa = sa.replaceAll("i", "l");
		sa = sa.replaceAll("v", "b");
		sa = sa.replaceAll("s", "p");
		sa = sa.replaceAll("c", "d");
		sa = sa.replaceAll("h", "g");
		sa = sa.replaceAll("f", "a");
		return sa;
	}

	/**
	 * e l d inside outside
	 * 
	 * @param sa
	 * @return
	 */
	public String getSAV3(String sa) {
		sa = sa.replaceAll("k", "e");
		sa = sa.replaceAll("i", "l");
		// sa = sa.replaceAll("v", "b");
		// sa = sa.replaceAll("s", "p");
		sa = sa.replaceAll("c", "d");
		// sa = sa.replaceAll("h", "g");
		// sa = sa.replaceAll("f", "a");
		return sa;
	}

	public String getSA(StrucState[] pr, String strAccs)
			throws StructureException {

		// StrucState[] pr = new StrucState[atoms.length];
		// for (int i = 0; i < atoms.length; i++) {
		//
		// pr[i] = new StrucState();
		// }

		StringBuilder sabBuilder = new StringBuilder();
		// for (int i = 0; i < atoms.length - 1; i++) {
		//
		// Group groupa = atoms[i].getGroup();
		// Group groupb = atoms[i + 1].getGroup();
		//
		// Atom a_N = (Atom) groupa.getAtomByPDBname(" N  ").clone();
		// Atom a_CA = (Atom) groupa.getAtomByPDBname(" CA ").clone();
		// Atom a_C = (Atom) groupa.getAtomByPDBname(" C  ").clone();
		// // Atom a_O = (Atom) groupa.getAtomByPDBname(" O  ").clone();
		// Atom b_N = (Atom) groupb.getAtomByPDBname(" N  ").clone();
		// Atom b_CA = (Atom) groupb.getAtomByPDBname(" CA ").clone();
		// Atom b_C = (Atom) groupb.getAtomByPDBname(" C  ").clone();
		// // Atom b_O = (Atom) groupb.getAtomByPDBname(" O  ").clone();
		//
		// double phi = Calc.torsionAngle(a_C, b_N, b_CA, b_C);
		//
		// double psi = Calc.torsionAngle(a_N, a_CA, a_C, b_N);
		//
		// double omega = Calc.torsionAngle(a_CA, a_C, b_N, b_CA);
		//
		// pr[i + 1].setPhi(phi);
		// pr[i].setPsi(psi);
		// pr[i].setOmega(omega);
		//
		// }

		for (int i = 0; i < pr.length; i++) {

			double phi = pr[i].getPhi();
			double psi = pr[i].getPsi();
			// double omega = pr[i].getOmega();
			char acc = '0';
			// if (i < strAccs.length())
			acc = strAccs.charAt(i);

			if (40.0 <= phi && phi <= 160.0 && -180.0 <= psi && psi <= -40.0) {
				if (acc == '0')// inside
					sabBuilder.append("k");
				else
					sabBuilder.append("e");
			}

			else if (40.0 <= phi && phi <= 160.0 && 145.0 <= psi
					&& psi <= 180.0) {
				if (acc == '0')// inside
					sabBuilder.append("k");
				else
					sabBuilder.append("e");
			} else if (30.0 <= phi && phi <= 100.0 && -20.0 <= psi
					&& psi <= 90.0) {
				if (acc == '0')// inside
					sabBuilder.append("i");
				else
					sabBuilder.append("l");
			} else if (-180.0 <= phi && phi <= -90.0 && 90.0 <= psi
					&& psi <= 180.0) {
				if (acc == '0')// inside
					sabBuilder.append("v");
				else
					sabBuilder.append("b");
			} else if (-180.0 <= phi && phi <= -90.0 && -180.0 <= psi
					&& psi <= -165.0) {
				if (acc == '0')// inside
					sabBuilder.append("v");
				else
					sabBuilder.append("b");
			} else if (160.0 <= phi && phi <= 180.0 && 90.0 <= psi
					&& psi <= 180.0) {
				if (acc == '0')// inside
					sabBuilder.append("v");
				else
					sabBuilder.append("b");
			} else if (160.0 <= phi && phi <= 180.0 && -180.0 <= psi
					&& psi <= -165.0) {
				if (acc == '0')// inside
					sabBuilder.append("v");
				else
					sabBuilder.append("b");
			} else if (-90.0 <= phi && phi <= -30.0 && 90.0 <= psi
					&& psi <= 180.0) {
				if (acc == '0')// inside
					sabBuilder.append("s");
				else
					sabBuilder.append("p");
			} else if (-90.0 <= phi && phi <= -30.0 && -180.0 <= psi
					&& psi <= -165.0) {
				if (acc == '0')// inside
					sabBuilder.append("s");
				else
					sabBuilder.append("p");
			} else if (180.0 <= phi && phi <= -50.0 && 20.0 <= psi
					&& psi <= 90.0) {
				if (acc == '0')// inside
					sabBuilder.append("c");
				else
					sabBuilder.append("d");
			} else if (-160.0 <= phi && phi <= -50.0 && -25.0 <= psi
					&& psi <= 20.0) {
				if (acc == '0')// inside
					sabBuilder.append("h");
				else
					sabBuilder.append("g");
			} else if (-140.0 <= phi && phi <= -20.0 && -80.0 <= psi
					&& psi <= -25.0) {
				if (acc == '0')// inside
					sabBuilder.append("f");
				else
					sabBuilder.append("a");
			} else {

				// if (-30 <= omega && omega <= 150.0)
				sabBuilder.append("x");
			}
		}

		char[] sab = sabBuilder.toString().toCharArray();

		for (int i = 0; i < sab.length - 1; i++) {

			if (-30 <= pr[i].getOmega() && pr[i].getOmega() <= 150.0) {

				if (sab[i] == 'x' && sab[i + 1] == 'x') {
					sab[i] = 'y';
					sab[i + 1] = 'y';

				}
			}

		}

		return new String(sab);

	}
}
