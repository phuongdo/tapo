package org.cnrs.crbm.lib.cm;

import java.io.BufferedWriter;
import java.io.FileWriter;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.utils.PdbTools;

public class ContactSignals {

	public static void main(String[] args) throws StructureException {

		String pdbCode = args[0];
		String pdbChain = args[1];
		// String pdbFile = PdbTools.downloadPDB(pdbCode);
		Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
		Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
				.getAtomCAArray(structure.getChainByPDB(pdbChain)));
		DSSP dssp = new DSSP(pdbCode);
		ContactMap contactMap = new ContactMap(cbAtoms, 7.0, dssp.filterSS(dssp
				.getSS(pdbChain, cbAtoms)));
		contactMap.setPdbCode(pdbCode + "_" + pdbChain);
		double[] histogram = contactMap.getHistogram();

		// write to matlab format???
		// C:\Program Files\matlab\toolbox\funts

		String matlab_dir = "C:\\Users\\CRBM\\workspace\\Matlab\\data\\";
		String outputfile = matlab_dir + pdbCode + "_" + pdbChain + ".sig";

		try {
			FileWriter fstream = new FileWriter(outputfile);
			BufferedWriter out = new BufferedWriter(fstream);
			for (int i = 0; i < histogram.length; i++) {
				Atom atom = cbAtoms[i];
				int resiNo = atom.getGroup().getResidueNumber().getSeqNum();
				out.write(resiNo + "\t" + histogram[i] + "\n");
			}

			// i love this current keyboard?

			out.close();
			fstream.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
