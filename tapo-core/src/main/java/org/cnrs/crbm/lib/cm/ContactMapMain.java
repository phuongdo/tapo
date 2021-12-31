package org.cnrs.crbm.lib.cm;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.Options;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.utils.PdbTools;

public class ContactMapMain {

	public static void main(String[] args) {

		final CommandLineParser cmdLineGnuParser = new GnuParser();
		final Options options = constructOptions();
		CommandLine commandLine;
		String pdbFile = null, pdbChain = null, proteinCode = null;
		double cutoff = 0.0;
		boolean isLoadDB = false;
		try {
			commandLine = cmdLineGnuParser.parse(options, args);

			// p or i
			if (commandLine.hasOption("p")) {
				proteinCode = commandLine.getOptionValue("p");
				isLoadDB = true;
			} else {
				pdbFile = commandLine.getOptionValue("i");
				isLoadDB = false;
			}
			if (commandLine.hasOption("c")) {
				pdbChain = commandLine.getOptionValue("c");

			}
			if (commandLine.hasOption("cutoff")) {
				cutoff = Double.parseDouble(commandLine
						.getOptionValue("cutoff"));

			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}

		ContactMapMain cmTest = new ContactMapMain();
		if (!isLoadDB)
			cmTest.run(pdbFile, pdbChain, cutoff);

		else {
			// get pdbFile from temporary directory

			pdbFile = PdbTools.downloadPDB(proteinCode);
			cmTest.run(pdbFile, pdbChain, cutoff);
		}
	}

	private static Options constructOptions() {
		final Options options = new Options();
		options.addOption("p", true, "PDB code").addOption("c", true, "Chain")
				.addOption("cutoff", true, "cutoff")
				.addOption("i", true, "input ex. input/1BPO");

		return options;
	}

	public void run(String pdbFile, String pdbChain, double cutoff) {

		// String pdb1 = "input/1BPO.pdb";

		PDBFileReader pdbreader = new PDBFileReader();
		String pdbCode = "";
		Atom[] atomSet1 = null;

		try {
			Structure struc1 = pdbreader.getStructure(pdbFile);
			pdbCode = struc1.getPDBCode();
			Chain chain = struc1.getChainByPDB(pdbChain);
			atomSet1 = PdbTools.getAtomCBArray(StructureTools
					.getAtomCAArray(chain));
			// System.out.println(atomSet1.length);

		} catch (Exception e) {
			e.printStackTrace();
		}
		DSSP dssp = new DSSP(pdbCode);
		ContactMapIO cmIO = new ContactMapIO();
		ContactMap contactMap = new ContactMap(atomSet1, cutoff, dssp.getSS(
				pdbChain, atomSet1));
		contactMap.setPdbCode(pdbCode + "_" + pdbChain);
		// cmIO.saveCMtoPNG(contactMap);
		// cmIO.saveCMtoFile(contactMap);
		cmIO.saveHistogramtoPNG(contactMap);

	}
}
