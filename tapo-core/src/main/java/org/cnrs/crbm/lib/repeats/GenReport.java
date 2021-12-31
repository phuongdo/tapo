package org.cnrs.crbm.lib.repeats;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class GenReport {

	public static void main(String[] args) {

		// type repeat or non-repeats
		GenReport genReport = new GenReport();
		String TYPE = args[0];
		List<String> listPros = genReport.getProtein(TYPE);
		genReport.exportLatexFile(listPros, TYPE);

	}

	public void exportLatexFile(List<String> strs, String type) {

		String matlab_dir = "C:\\Users\\CRBM\\workspace\\Matlab\\";
		String outputfile = matlab_dir + type + "_reports.tex";

		StringBuilder builder = new StringBuilder();
		builder.append("\\documentclass[10pt,a4paper]{article}\n");
		builder.append("\\usepackage[margin=1cm]{geometry}\n");
		builder.append("\\usepackage{graphicx}\n");
		builder.append("\\usepackage{subfig}\n");
		builder.append("\\begin{document}\n");

		int i = 0;

		for (String str : strs) {
			// get db code and chain
			String pdbCode = str.substring(0, 4).toUpperCase();
			String pdbChain = str.substring(4, 5);
			// System.out.println(pdbCode + " " + pdbChain);
			if (i % 3 == 0)
				builder.append("\\clearpage\n");
			builder.append("\\begin{figure}\n");
			String imgSignal = "output/" + type + "/vect/" + pdbCode + "_"
					+ pdbChain + "_vect.png";
			String imgFourier = "output/" + type + "/vect/" + pdbCode + "_"
					+ pdbChain + "_fourier.png";

			builder.append("\\subfloat[RMSD signals]{{\\includegraphics[width=9cm]{"
					+ imgSignal + "} }}\n");
			builder.append("\\qquad\n");
			builder.append("\\subfloat[Fourier Tranform]{{\\includegraphics[width=9cm]{"
					+ imgFourier + "} }}\n");
			builder.append("\\caption{Protein " + pdbCode + " " + pdbChain
					+ "}");
			builder.append("\\label{fig:example}\n");
			builder.append("\\end{figure}\n");
			i++;
		}

		builder.append("\\end{document}\n");

		try {
			FileWriter fstream = new FileWriter(outputfile);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(builder.toString());
			out.close();
			fstream.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public List<String> getProtein(String type) {

		List<String> strs = new ArrayList<String>();
		BufferedReader br = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader("input/" + type));
			while ((sCurrentLine = br.readLine()) != null) {
				strs.add(sCurrentLine);

			}

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return strs;
	}

}
