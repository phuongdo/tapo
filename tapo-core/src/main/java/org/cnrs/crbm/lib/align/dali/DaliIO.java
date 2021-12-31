package org.cnrs.crbm.lib.align.dali;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import org.biojava.nbio.structure.*;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.utils.AlignUtils;

public class DaliIO {

	public AlignChainDali readOutput(String dir) {

		AlignChainDali alignChain = new AlignChainDali();
		try {

			// Open the file that is the first
			// command line parameter
			FileInputStream fstream = new FileInputStream(dir);
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;

			// read header
			strLine = br.readLine();
			strLine = br.readLine();

			// read meta pairAlign

			while ((strLine = br.readLine()) != null && !strLine.equals("")) {
				// Print the content on the console
				// System.out.println(strLine);

				if (strLine.contains("#")) {

					break;
				}

				else {
					// just read one
					int index = Integer
							.parseInt(strLine.substring(0, 4).trim());
					if (index == 1) {

						// System.out.println(strLine);

						double rmsd = Double.parseDouble(strLine.substring(20,
								23).trim());

						int algnLengh = Integer.parseInt(strLine.substring(26,
								28).trim());
						// System.out.println(rmsd);
						alignChain.setRmsd(rmsd);
						alignChain.setChainLength(algnLengh);
					}

				}

			}

			List<Integer> pair1 = new ArrayList<Integer>();
			List<Integer> pair2 = new ArrayList<Integer>();

			while ((strLine = br.readLine()) != null && !strLine.equals("")) {
				// Print the content on the console
				// if (strLine.equals(""))
				// continue;

				if (!strLine.contains("#")) {

					int index = Integer
							.parseInt(strLine.substring(0, 4).trim());

					int lastEnd = 0;
					if (index == 1) {

						// System.out.println(strLine);
						int pair1Start = Integer.parseInt(strLine.substring(20,
								25).trim());
						int pair1End = Integer.parseInt(strLine.substring(27,
								31).trim());

						int pair2Start = Integer.parseInt(strLine.substring(36,
								40).trim());
						int pair2End = Integer.parseInt(strLine.substring(42,
								47).trim());
						// System.out.println(pair1Start + "-" + pair1End +
						// " = "
						// + pair2Start + "-" + pair2End);
						for (int i = 0; i < pair1End - pair1Start + 1; i++) {

							alignChain.getAlignPairs().put(pair1Start + i - 1,
									pair2Start + i - 1);
							pair1.add(pair1Start + i - 1);
							pair2.add(pair2Start + i - 1);

						}

					}

				}

			}

			int[][] multi = new int[2][pair1.size()];
			for (int i = 0; i < pair1.size(); i++) {
				multi[0][i] = pair1.get(i);
				multi[1][i] = pair2.get(i);
			}

			int gaps = 0;
			List<List<Integer>> multiples = AlignUtils.getMultipleAligns(multi,
					2, pair1.size());

			for (List<Integer> alist : multiples) {
				for (Integer k : alist) {
					if (k == -1)
						gaps++;
				}
			}
			alignChain.setGaps(gaps);
			// Close the input stream
			in.close();
		} catch (Exception e) {// Catch exception if any
			e.printStackTrace();
		}
		return alignChain;

	}

	/**
	 * Write a pdb file from structure object.
	 * 

	 * @return
	 */

	public Map<Integer, Integer> writePDB(Atom[] setAtom, String fileName) {

		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		try {
			// create new Structure
			Structure struc = new StructureImpl();
			Chain c = new ChainImpl();

			for (int i = 0; i < setAtom.length; i++) {
				Atom atom = setAtom[i];
				map.put(atom.getGroup().getResidueNumber().getSeqNum(), i);
				c.addGroup(atom.getGroup());
			}
			struc.addChain(c);

			// write to pdb file
			String outputfile = Dir.DALI_TMP_DIR + "/" + fileName + ".pdb";
			FileOutputStream out = new FileOutputStream(outputfile);
			PrintStream p = new PrintStream(out);
			p.println(struc.toPDB());
			p.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

		return map;
	}
}
