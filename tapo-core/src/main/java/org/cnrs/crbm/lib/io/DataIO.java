package org.cnrs.crbm.lib.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.cnrs.crbm.lib.cdhit.CDHIT;

public class DataIO {

	public static void writeToFile(String content, String outputDir) {

		try {
			// Create file
			FileWriter fstream = new FileWriter(outputDir);
			BufferedWriter out = new BufferedWriter(fstream);
			// System.out.println(content);
			out.write(content);
			out.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

	}

	public static String readFile(String fileDir) {
		StringBuilder builder = new StringBuilder();
		try {

			// Open the file that is the first
			// command line parameter
			FileInputStream fstream = new FileInputStream(fileDir);
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			while ((strLine = br.readLine()) != null) {
				builder.append(strLine + "\n");

			}
			// Close the input stream
			in.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

		return builder.toString();
	}

	public static List<String> readLines(String fileDir) {
		List<String> lines = new ArrayList<String>();
		try {

			// Open the file that is the first
			// command line parameter
			FileInputStream fstream = new FileInputStream(fileDir);
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			while ((strLine = br.readLine()) != null) {

				if(strLine.startsWith("#"))
					continue;
				if(strLine.trim().equals(""))
					continue;
				lines.add(strLine);

			}
			// Close the input stream
			in.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

		return lines;
	}

	public static List<Row> getListProteinFromFile(String fileDir) {

		List<Row> rows = new ArrayList<Row>();

		List<String> lines = readLines(fileDir);
		Set<String> pdbs = new HashSet<String>();
		for (String line : lines) {

			if(line.startsWith("#"))
				continue;
			pdbs.add(line.substring(0, 6));

		}

		for (String pdb : pdbs) {
			Row row = new Row();
			row.setProtein(pdb);
			rows.add(row);
		}

		return rows;

	}

	public static List<Row> getListProteinFromFilePdb(String fileDir) {
		DecimalFormat df = new DecimalFormat("0.00");
		List<Row> rows = new ArrayList<Row>();

		List<String> lines = readLines(fileDir);
		int size = lines.size();
		int i = 0;
		// System.out.println("loading");
		for (String line : lines) {
			String pdb = line.toLowerCase();

			double percent = (double) i * 100 / size;
			System.out.print("\r loading " + df.format(percent) + " %");
			CDHIT cdhit = new CDHIT(pdb);

			List<String> chains = cdhit.getUniqueChain();

			for (String chain : chains) {
				Row row = new Row();
				row.setProtein(pdb + "_" + chain);
				// row.setProtein(line.substring(0, ));
				rows.add(row);
			}
			i++;
		}

		return rows;

	}
}
