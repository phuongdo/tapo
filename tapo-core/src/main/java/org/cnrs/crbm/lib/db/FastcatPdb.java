package org.cnrs.crbm.lib.db;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.zip.GZIPInputStream;

public class FastcatPdb {

	public static void main(String[] args) throws IOException {

		new FastcatPdb().extractFatcatPDB();
	}

	
	
	private void extractRepeatDB(){
		
		
	}
	
	private void extractFatcatPDB() throws IOException {

		String fileDir = "F:\\Storage\\fatcat_pdb\\";

		String gizFile = fileDir + "\\fatcat_rigid_pdb_all.txt.gz";
		String outFile = fileDir + "\\fatcat_rigid_pdb_all.csv";

		PrintWriter writer = new PrintWriter(outFile);
		// read and write data
		FileInputStream fin = new FileInputStream(gizFile);
		GZIPInputStream gzis = new GZIPInputStream(fin);
		InputStreamReader xover = new InputStreamReader(gzis);
		BufferedReader is = new BufferedReader(xover);
		String line;
		while ((line = is.readLine()) != null) {
			if (line.startsWith("#"))
				continue;
			String[] rows = line.split("\t");

			String name1 = rows[0];
			String name2 = rows[1];

			if (name1.startsWith("PDP") && name2.startsWith("PDP")) {
				// int leng1 = Integer.parseInt(rows[2]);
				// int leng2 = Integer.parseInt(rows[3]);
				// double pValue = Double.parseDouble(rows[4]);
				// double rmsd = Double.parseDouble(rows[5]);
				// double score = Double.parseDouble(rows[6]);
				// double pid = Double.parseDouble(rows[7]);
				// double tmScore = -1;
				// if (!rows[8].equals("nu"))
				// tmScore = Double.parseDouble(rows[8]);
				// int cov1 = Integer.parseInt(rows[9]);
				// int cov2 = Integer.parseInt(rows[10]);

				name1 = name1.substring(4, 9);
				name2 = name2.substring(4, 9);

				writer.write(name1 + ";" + name2 + ";" + rows[2] + ";"
						+ rows[3] + ";" + rows[4] + ";" + rows[5] + ";"
						+ rows[6] + ";" + rows[7] + ";" + rows[8] + ";"
						+ rows[9] + ";" + rows[10] + "\n");
			}

		}

		writer.close();
		is.close();
		xover.close();
		gzis.close();
		fin.close();

	}
}
