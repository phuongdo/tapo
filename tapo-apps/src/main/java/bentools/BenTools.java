package bentools;

import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.raphael.Raphael;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.util.List;

/**
 * Created by pdoviet on 9/8/2015.
 */
public class BenTools {

    public static void main(String[] args) throws Exception {

        new BenTools().benchmark();

    }

    public void benchmark() throws Exception {
//        List<String> datasetTRs = DataIO.readLines("data/benchmarkset/pub/Dataset-TR.txt");
//        for (String row : datasetTRs) {
//            try {
//                processRaphael(row, "1");
//            } catch (Exception ex) {
//                ex.printStackTrace();
//            }
//        }
        List<String> datasetNonTRs = DataIO.readLines("data/benchmarkset/pub/Dataset-non-TR-v.1.1.0.txt");
        for (String row : datasetNonTRs) {
            try {
                processRaphael(row, "0");
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

    }

    public void processRaphael(String row, String label) throws Exception {
//        System.out.println(row);

        String[] data = row.split("\t");
        String pdb = data[0];
        String pdbCode = pdb.substring(0, 4);
        String pdbChain = pdb.substring(5, 6);
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        Atom[] atoms = repeatFinder.getAtoms();
        Raphael raphael = new Raphael(atoms);
        System.out.println(row + "\t" + NumberFormatUtils.format(raphael.getTotalScore()) + "\t" + label);


    }

}
