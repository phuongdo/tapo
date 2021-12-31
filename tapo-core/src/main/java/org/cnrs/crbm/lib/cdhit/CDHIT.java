package org.cnrs.crbm.lib.cdhit;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import org.biojava.nbio.structure.Structure;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.utils.ExecUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

public class CDHIT {

	String cdHitOutFile;

	String proCode;
	double c;
	int n;

	/**
	 * get from cache
	 */
	public CDHIT(String pdbCode) {
		this.c = 0.9;
		this.n = 5;
		this.proCode = pdbCode;
		this.cdHitOutFile = Dir.FASTA_LOCAL + "/" + proCode.substring(1, 3)
				+ "/" + proCode + ".cdhit";

		try {

			File f = new File(this.cdHitOutFile);
			if (!f.exists()) {

				Structure structure = PdbTools
						.getStructureFromLocalPdb(proCode);
				String pdbFileFasta = PdbTools.insertToFastaDB(structure);
				calBashshell(pdbFileFasta);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * genereate cdhit cluster of similarity sequence of one protein.
	 * 
	 * @param structure
	 * @param c
	 * @param n
	 */
	public CDHIT(Structure structure, double c, int n) {

		this.c = c;
		this.n = n;
		this.proCode = structure.getPDBCode().toLowerCase();
		this.cdHitOutFile = Dir.FASTA_LOCAL + "/" + proCode.substring(1, 3)
				+ "/" + proCode + ".cdhit";

		try {

			File f = new File(this.cdHitOutFile);
			if (!f.exists()) {
				String pdbFileFasta = PdbTools.insertToFastaDB(structure);
				calBashshell(pdbFileFasta);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	// simple
	// public CDHIT(String proCode, double c, int n) {
	// this.proCode = proCode;
	// this.cdHitOutFile = Dir.CDHIT_TMP_DIR + "/" + proCode + ".cdhit";
	// this.c = c;
	// this.n = n;
	// String pdbFileFasta = PdbTools.downloadFastaFromServer(proCode);
	// try {
	//
	// File f = new File(this.cdHitOutFile);
	// if (!f.exists()) {
	// calBashshell(pdbFileFasta);
	// }
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	// }

	private void calBashshell(String pdbFileFasta) {

		String shellCommand = Dir.CDHIT_EXECUTABLE + " -i " + pdbFileFasta
				+ " -o " + this.cdHitOutFile + " -c " + c + " -n " + n;
		// System.out.println(shellCommand);
		String strCmdOut = ExecUtils.execShellCmdLinux(shellCommand);

	}

	public List<String> getUniqueChain() {

		List<String> uniqueChains = new ArrayList<String>();
		try {
			// Open the file that is the first
			// command line parameter
			FileInputStream fstream = new FileInputStream(this.cdHitOutFile);
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String line;

			// Read File Line By Line

			while ((line = br.readLine()) != null) {
				Pattern p = Pattern.compile(".*\\:(.*)\\|.*\\|.*\\|.*");
				Matcher m = p.matcher(line);
				if (m.find()) {
					// System.out.println(m.group(1));
					uniqueChains.add(m.group(1).trim());
				}
				// if(line.startsWith(">"))
				// System.out.println(line);
			}

			// Close the input stream
			in.close();
		} catch (Exception e) {// Catch exception if any
			// System.err.println("Error: " + e.getMessage());
		}
		return uniqueChains;
	}

	public void clearTemporaryFiles() {

		new File(this.cdHitOutFile).delete();
	}

	public static void main(String[] args) {

		String pdbCode = "1WCS";
		if (args.length > 0)
			pdbCode = args[0];

		pdbCode = pdbCode.toLowerCase();
		Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
		CDHIT cdhit = new CDHIT(structure, 0.9, 5);
		List<String> uniqueChains = cdhit.getUniqueChain();
		for (String chain : uniqueChains) {
			System.out.println(chain);
		}

		// cdhit.clearTemporaryFiles();

	}
}
