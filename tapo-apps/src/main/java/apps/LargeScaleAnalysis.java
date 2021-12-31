package apps;

import java.io.*;
import java.util.List;
import java.util.Set;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.repeats.CombineScore;
import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.sadb.Sequence3D;
import org.cnrs.crbm.lib.utils.ProgressBar;

public class LargeScaleAnalysis {
    public static void main(String[] args) throws Exception {

        StringBuilder builder = new StringBuilder();
        List<String> lines = DataIO.readLines("data/resultLargeScale/PDB_26June2015.o");
        for (String line : lines) {
            try {
                if (line.contains(">")) {
                    //line = line.substring(1,line.length());
                    String pdb = line.split("\\|")[0].replace(">", "");
                    String datas = line.split("\\|")[1];
                    CombineScore combineScore = new CombineScore(datas);
                    builder.append(pdb + "," + combineScore.toCSVFormat() + "\n");


                }
            } catch (Exception ex) {
                ex.printStackTrace();
                //System.out.println(line);
            }


        }

        DataIO.writeToFile(builder.toString(), "output/PDB_26June2015.out");


//        ReadFasta readFasta = new ReadFasta();
//
//
//        // read pdbs from none redundancy db
//        Set<String> pdbs = readFasta.getAllPdbsWithChain("data/PDB_July_1_2011.txt");
//
//        for (String pdb : pdbs) {
//
//            String pdbCode = pdb.substring(0, 4).toLowerCase();
//            String pdbChain = pdb.substring(5, 6);
//            System.out.println(pdbCode + "_" + pdbChain);
//
//
//        }

//		try {
//
//
//
//
//			// Open the file that is the first
//			// command line parameter
//			// String input = "input/RepeatData140618.csv";
//			// String output = "output/francois_reports.html";
//			// String input = "input/repeat-raf";
//			// String output = "output/raf-db-reports.html";
//			String input = args[0];
//			String output = args[1];
//			ProteinCSVReader csvReader = new ProteinCSVReader();
//			List<Row> rows = csvReader.getData(input);
//			StringBuilder builder = new StringBuilder();
//
//			builder.append("<table>");
//			builder.append("<tr><td>protein</td><td>Tag</td><td>Conclusion</td><td>RMSD</td><td>T-REKs 3D</td><td>TRUST 3D</td><td>TRUST 3D(*)</td><td>VECTORS</td><td>VECTORS-2ndS</td><td>CONTACT MAP</td><td>RAPHEAL</td><tr/>");
//
//			for (Row row : rows) {
//				// System.out.println(strLine);
//
//				String pdbCode = row.getPdbCode();
//				String pdbChain = row.getPdbChain();
//
//				System.out.println(pdbCode + "_" + pdbChain);
//
//				try {
//					RepeatFinder repeatFinder = new RepeatFinder(pdbCode,
//							pdbChain);
//					repeatFinder.setMode(FinderMode.CONSOLE);
//					repeatFinder.findRepeat();
//					// repeatFinder.debug();
//					builder.append(repeatFinder
//							.getCompareResultsBetweenEachMethos(row) + "\n");
//
//				} catch (Exception ex) {
//					// builder.append("<tr><td></td><tr>");
//					ex.printStackTrace();
//
//				}
//
//			}
//
//			builder.append("</table>");
//			DataIO.writeToFile(builder.toString(), output);
//
//		} catch (Exception e) {// Catch exception if any
//			System.err.println("Error: " + e.getMessage());
//		}
    }


}
