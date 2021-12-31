package apps;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.util.List;

/**
 * Created by pdoviet on 7/2/2015.
 */
public class UpdatePDB {

    public static void main(String[] args) {

        String fileOutDir = "data/pdbList.in";
        if (args.length > 0) {
            fileOutDir = args[0];

        }
        try {

            System.out.println("Starting read : " + fileOutDir);
            StringBuffer buffer = new StringBuffer();
            List<String> pdbs = DataIO.readLines(fileOutDir);
            System.out.println("Load job input with pdb size : " + pdbs.size());
            System.out.println("Download data...");
            ProgressBar bar = new ProgressBar();
            int process = 1;
            int sizeOfProcess = pdbs.size();
            RepeatFinder repeatFinder = null;
            //ExperRMSD finder = null;
            for (String pdb : pdbs) {
                // process and get output
                try {
                    String pdbCode = pdb.substring(0, 4);
                    String pdbChain = pdb.substring(5, 6);
                    repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                    repeatFinder.setMode(FinderMode.CONSOLE);
                    repeatFinder.setNrOfProcessors(1);
                    bar.update(process, sizeOfProcess);
                    process++;

                } catch (Exception e) {

                }
            }

            System.out.println("Done!!!");

        } catch (Exception e) {

        }

    }
}
