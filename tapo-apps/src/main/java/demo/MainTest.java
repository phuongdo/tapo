package demo;

import com.apporiented.algorithm.clustering.AverageLinkageStrategy;
import com.apporiented.algorithm.clustering.Cluster;
import com.apporiented.algorithm.clustering.ClusteringAlgorithm;
import com.apporiented.algorithm.clustering.DefaultClusteringAlgorithm;
import nptr.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.biojava.nbio.structure.*;
import org.cnrs.crbm.lib.utils.PdbTools;

import javax.swing.*;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class MainTest implements Runnable {

	public static void main(String[] args) throws Exception {

		// test for logging

		String[] names = new String[] { "O1", "O2", "O3", "O4", "O5", "O6" };
		double[][] distances = new double[][] { { 0, 1, 9, 7, 11, 14 },
				{ 1, 0, 4, 3, 8, 10 }, { 9, 4, 0, 9, 2, 8 },
				{ 7, 3, 9, 0, 6, 13 }, { 11, 8, 2, 6, 0, 10 },
				{ 14, 10, 8, 13, 10, 0 } };

		ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
		Cluster cluster = alg.performClustering(distances, names,
				new AverageLinkageStrategy());

		System.out.println(cluster.getChildren().get(1).getChildren().get(0)
				.getTotalDistance());

		//
		// DendrogramPanel dp = new DendrogramPanel();
		// dp.setModel(cluster);
		//
		// // 1. Create the frame.
		// JFrame frame = new JFrame("FrameDemo");
		// frame.setSize(300, 400);
		//
		// // 2. Optional: What happens when the frame closes?
		// frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		//
		// // 3. Create components and put them in the frame.
		// // ...create emptyLabel...
		// frame.getContentPane().add(dp, BorderLayout.CENTER);
		//
		// // 4. Size the frame.
		// frame.pack();
		//
		// // 5. Show it.
		// frame.setVisible(true);

		// MainTest test = new MainTest();
		// ProteinCSVReader csvReader = new ProteinCSVReader();
		// List<Row> rows = csvReader.getData("input/RepeatData140723.csv");
		// System.out.println(rows.size());
		//
		// for (Row row : rows) {
		// // System.out.println(strLine);
		//
		// String pdbCode = row.getPdbCode();
		// String pdbChain = row.getPdbChain();
		//
		// System.out.println(pdbCode + "_" + pdbChain);
		// }
		// // test.test3();
		// test.test1();
		// test.test4();
		// test.test4();
		// String[] argurments = new String[] {
		// "-infile=C:\\Users\\CRBM\\Desktop\\2F8X.fasta.txt" };
		// NPTR nptr = new NPTR(argurments);

	}

	public void testLogging() {

	}

	public void test5() {

		try {
			String input = "input/repeat-raf";
			String output = "output/raf-db-reports.html";
			FileInputStream fstream = new FileInputStream(input);
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;

			while ((strLine = br.readLine()) != null) {
				// System.out.println(strLine);

				String pdbCode = strLine.substring(0, 4).toUpperCase();
				String pdbChain = strLine.substring(4, 5).toUpperCase();
				if (!pdbChain.equals("0"))
					System.out.println(pdbCode + "_" + pdbChain);
				else
					System.out.println(pdbCode + "_A");

			}
			in.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}

	public void test4() {

		char[] atr = "EKLIBVPSDCGHAFXYOWMN".toCharArray();

		// print first
		System.out.print("   ");
		for (int i = 0; i < atr.length; i++) {
			System.out.print(atr[i] + "  ");
		}

		System.out.println();
		for (int i = 0; i < atr.length; i++) {

			System.out.print(atr[i] + " ");

			for (int j = 0; j < atr.length; j++) {

				if (atr[i] == atr[j] && atr[i] != 'X' && atr[i] != 'Y') {
					if (!"OWMN".contains(atr[i] + ""))
						System.out.print(" 5" + " ");
					else
						System.out.print(" 10" + " ");
				}
				// e(k)
				else if ((atr[i] == 'E' && atr[j] == 'K')
						|| (atr[i] == 'K' && atr[j] == 'E'))
					System.out.print(" 4" + " ");

				// l(i)
				else if ((atr[i] == 'L' && atr[j] == 'I')
						|| (atr[i] == 'I' && atr[j] == 'L'))
					System.out.print(" 4" + " ");

				// b(v) ->p(s)

				else if ((atr[i] == 'B' && atr[j] == 'V')
						|| (atr[i] == 'V' && atr[j] == 'B'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'P' && atr[j] == 'S')
						|| (atr[i] == 'S' && atr[j] == 'P'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'B' && atr[j] == 'P')
						|| (atr[i] == 'P' && atr[j] == 'B'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'V' && atr[j] == 'S')
						|| (atr[i] == 'S' && atr[j] == 'V'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'B' && atr[j] == 'S')
						|| (atr[i] == 'S' && atr[j] == 'B'))
					System.out.print(" 2" + " ");

				else if ((atr[i] == 'P' && atr[j] == 'V')
						|| (atr[i] == 'V' && atr[j] == 'P'))
					System.out.print(" 2" + " ");

				// a(f)->g(h)

				else if ((atr[i] == 'A' && atr[j] == 'G')
						|| (atr[i] == 'G' && atr[j] == 'A'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'A' && atr[j] == 'F')
						|| (atr[i] == 'F' && atr[j] == 'A'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'G' && atr[j] == 'H')
						|| (atr[i] == 'H' && atr[j] == 'G'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'F' && atr[j] == 'H')
						|| (atr[i] == 'H' && atr[j] == 'F'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'A' && atr[j] == 'H')
						|| (atr[i] == 'H' && atr[j] == 'A'))
					System.out.print(" 2" + " ");

				else if ((atr[i] == 'F' && atr[j] == 'G')
						|| (atr[i] == 'G' && atr[j] == 'F'))
					System.out.print(" 2" + " ");

				// g(h)->d(c)

				else if ((atr[i] == 'G' && atr[j] == 'D')
						|| (atr[i] == 'D' && atr[j] == 'G'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'H' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'H'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'G' && atr[j] == 'H')
						|| (atr[i] == 'H' && atr[j] == 'G'))
					System.out.print(" 4" + " ");
				else if ((atr[i] == 'D' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'D'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'G' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'G'))
					System.out.print(" 2" + " ");
				else if ((atr[i] == 'H' && atr[j] == 'D')
						|| (atr[i] == 'D' && atr[j] == 'H'))
					System.out.print(" 2" + " ");

				// b(v)->d(c)
				else if ((atr[i] == 'B' && atr[j] == 'D')
						|| (atr[i] == 'D' && atr[j] == 'B'))
					System.out.print(" 2" + " ");

				else if ((atr[i] == 'V' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'V'))
					System.out.print(" 2" + " ");

				else if ((atr[i] == 'B' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'B'))
					System.out.print(" 1" + " ");
				else if ((atr[i] == 'D' && atr[j] == 'V')
						|| (atr[i] == 'V' && atr[j] == 'D'))
					System.out.print(" 1" + " ");

				else if ((atr[i] == 'B' && atr[j] == 'V')
						|| (atr[i] == 'V' && atr[j] == 'B'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'D' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'D'))
					System.out.print(" 4" + " ");

				// p(s)->d(c)
				else if ((atr[i] == 'P' && atr[j] == 'D')
						|| (atr[i] == 'D' && atr[j] == 'P'))
					System.out.print(" 2" + " ");

				else if ((atr[i] == 'S' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'S'))
					System.out.print(" 2" + " ");

				else if ((atr[i] == 'P' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'P'))
					System.out.print(" 1" + " ");
				else if ((atr[i] == 'D' && atr[j] == 'S')
						|| (atr[i] == 'S' && atr[j] == 'D'))
					System.out.print(" 1" + " ");

				else if ((atr[i] == 'P' && atr[j] == 'S')
						|| (atr[i] == 'S' && atr[j] == 'P'))
					System.out.print(" 4" + " ");

				else if ((atr[i] == 'D' && atr[j] == 'C')
						|| (atr[i] == 'C' && atr[j] == 'D'))
					System.out.print(" 4" + " ");

				// o -n
				else if ((atr[i] == 'O' && atr[j] == 'N')
						|| (atr[i] == 'N' && atr[j] == 'O'))
					System.out.print(" 8" + " ");

				// b - v
				else if ((atr[i] == 'B' && atr[j] == 'V')
						|| (atr[i] == 'V' && atr[j] == 'B'))
					System.out.print(" 8" + " ");

				else {

					if ("OWMN".contains(atr[i] + "")
							|| "OWMN".contains(atr[j] + ""))
						System.out.print("-3" + " ");
					else
						System.out.print("-2" + " ");
				}

			}

			System.out.println();

		}

	}

	public void test3() throws Exception {
		// Parameters.setParamDefault("20");
		ParamField defaultParam = new ParamField();
		JFormattedTextField tLength = new JFormattedTextField();
		tLength.setText("Default");
		JFormattedTextField tPercent = new JFormattedTextField();
		tPercent.setText("20");
		Parameters.aParam.add(defaultParam);
		// Parameters.aParam.add(defaultParam);

		Parameters.threshold = 0.7;
		SeqRepeat seqRep = new SeqRepeat(
				"1satA",
				"HmlnlHllmBnnnnnHBlnBnllBlnnBlmHnlmnlllnllnllHnmBnlllBnBlBlnmBmBlBlBmBlBmBlmBnBHnnmnBBmmBnm");
		TRExplorer Trex = new TRExplorer(seqRep);
		OverlapManager aCopies = new OverlapManager();
		Trex.run();
		aCopies = Trex.getACopies();

		if (aCopies.size() >= 1) {

			for (int i = 0; i < aCopies.size(); i++) {

				if (((Repeat) aCopies.get(i)).isReal()) {
					System.out.println(((Repeat) aCopies.get(i)).toString()
							+ "\n");

					Repeat repeat = (Repeat) aCopies.get(i);
					System.out.println(repeat.getBeginPosition() + "_"
							+ repeat.getEndPosition());

				}
			}
			seqRep.clear();
			aCopies.clear();
		} else {
			System.out.println("repeat not found in sequence "
					+ seqRep.getDesc());
		}
	}



	public void test1() throws StructureException {

		Structure structure = PdbTools.getStructureFromFile("input/FAKE_A.pdb");
		// Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
		// .getAtomCAArray(structure.getChainByPDB(pdbChain)));
		//
		Atom[] atoms = StructureTools.getAtomCAArray(structure
				.getChainByPDB("C"));

		for (int i = 0; i < atoms.length - 4; i = i + 4) {

			System.out.println(i);
			Vector3D v1 = new Vector3D(atoms[i].getCoords());
			Vector3D v2 = new Vector3D(atoms[i + 1].getCoords());
			Vector3D v3 = new Vector3D(atoms[i + 2].getCoords());
			Vector3D v4 = new Vector3D(atoms[i + 3].getCoords());
			Vector3D v12 = v2.subtract(v1);
			Vector3D v34 = v4.subtract(v3);

			// System.out.println(Math.toDegrees(Vector3D.angle(v12, v34)));
			// Line line = new Line(v12, v34);
			System.out.println(Vector3D.distance(v12, v34));

		}

	}

	public void test2() throws Exception {

		String pdbCode = "3TLZ";
		PdbTools.downloadPDB("3TLZ");
		Structure structure = PdbTools.getStructureFromFile("input/3TLZ.pdb");
		Chain c = structure.getChainByPDB("A");

		// only the observed residues
		System.out.println(c.getAtomSequence());

		// print biological sequence
		System.out.println(c.getSeqResSequence());
	}

	// @Override
	public void run() {
		// TODO Auto-generated method stub

	}
}
