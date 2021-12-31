package org.cnrs.crbm.lib.annotation;

import java.util.List;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;

/**
 * Testing function for ProteinAnnotate class
 * 
 * @author pdoviet
 *
 */
public class ProteinAnnotateTest {

	public static void main(String[] args) throws Exception {

		// ProteinCSVReader csvReader = new ProteinCSVReader();
		// List<Row> rows = csvReader.getData("input/RepeatData140618.csv");
		// testGlobularProteins(rows);
		RepeatFinder repeatFinder = new RepeatFinder("4IRT", "A");
		repeatFinder.setMode(FinderMode.CONSOLE);
		ProteinAnnotate proteinAnnotate = new ProteinAnnotate();
		//
		// System.out.println(proteinAnnotate
		// .annotateDispersedProtein(repeatFinder.getAtoms()));

	}

	public static void testGlobularProteins(List<Row> rows) {
		ProteinAnnotate proteinAnnotate = new ProteinAnnotate();
		StringBuilder builder = new StringBuilder();
		for (Row row : rows) {
			// extracting features :)

			try {
				String pdbCode = row.getPdbCode();
				String pdbChain = row.getPdbChain();

				RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
				repeatFinder.setMode(FinderMode.CONSOLE);
				if (repeatFinder.getCmHistos().length == 0)
					continue;

				builder.append(pdbCode + "\t");
				// if (row.getByEyeTot() == 1)
				// builder.append("1");
				// else
				// builder.append("0");

				builder.append(row.getByEyeTot());

				System.out.print(pdbCode + "\t" + row.getByEyeTot() + "\t");

				boolean isGlobular = proteinAnnotate.annotateGlobularProtein(
						repeatFinder.getCmHistos(), 0, 0);

				builder.append("\t");
				if (isGlobular)
					builder.append("1");
				else
					builder.append("0");

				builder.append("\n");
			} catch (Exception ex) {

			}

		}

		DataIO.writeToFile(builder.toString(),
				"C:/Users/pdoviet/Desktop/report.txt");

	}

	public static void testLongHelix(List<Row> rows) {
		ProteinAnnotate proteinAnnotate = new ProteinAnnotate();
		StringBuilder builder = new StringBuilder();
		for (Row row : rows) {
			// extracting features :)

			try {
				String pdbCode = row.getPdbCode();
				String pdbChain = row.getPdbChain();

				RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
				repeatFinder.setMode(FinderMode.CONSOLE);
				if (repeatFinder.getStrSS().length() == 0)
					continue;

				if (row.getByEyeTot() == 4)
					builder.append("1");
				else
					builder.append("0");

				boolean isLongHelix = proteinAnnotate
						.annotateProteinWithLongHelix(repeatFinder.getStrSS(),
								0, 0);

				if (row.getByEyeTot() == 4 && !isLongHelix) {
					System.out.println(pdbCode);

					System.out.println(repeatFinder.getStrSS());
				}

				builder.append("\t");
				if (isLongHelix)
					builder.append("1");
				else
					builder.append("0");

				builder.append("\n");
			} catch (Exception ex) {

			}

		}

		// DataIO.writeToFile(builder.toString(),
		// "C:/Users/pdoviet/Desktop/report.txt");

	}
}
