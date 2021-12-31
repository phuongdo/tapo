package demo;

import org.cnrs.crbm.lib.io.DataIO;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by pdoviet on 8/19/2015.
 */
public class MetaFinder {

    public static void main(String[] args) {


        List<String> lines = DataIO.readLines("data/benchmarkset/pub/Dataset-TR.txt");

        Set<String> TRSet = new HashSet<String>();
        for (String line : lines) {
            TRSet.add(line);

        }
        lines = DataIO.readLines("data/metafinder/datasetWithLengthWithHeader");

        for (String line : lines) {

            String pdb = line.substring(0, 6);

            if (TRSet.contains(pdb) && line.contains("HHREPID")) {
                System.out.println(line);
            }
        }

        lines = DataIO.readLines("data/metafinder/repeatsDBcompleteDatasetWithLengthWithHeader");

        for (String line : lines) {
            String pdb = line.substring(0, 6);
            if (TRSet.contains(pdb) && line.contains("HHREPID")) {
                System.out.println(line);
            }
        }

    }
}
