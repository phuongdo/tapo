package apps;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.List;

public class ClusterJob {

    static Logger logger = LoggerFactory.getLogger(ClusterJob.class);

    public static void main(String[] args) {

        int jobid = 0;
        int noCore = 1;

        if (args.length > 0) {
            jobid = Integer.parseInt(args[0]);
            noCore = Integer.parseInt(args[1]);
        }
        String fileOutDir = Dir.CLUSTER_WORKING_DIR + "/cls_output/output."
                + jobid;


        boolean isExists = new File(fileOutDir).exists();

        // System.out.println(fileOutDir);
        if (isExists) {

            System.out.println("job output " + jobid + " exists!!!");
            System.exit(1);
        }

        try {

//            logger.info("Starting job : " + jobid);
            // System.out.println(Dir.FASTA_LOCAL + "/job_array/input" + jobid);
            StringBuffer buffer = new StringBuffer();
            List<String> pdbs = DataIO.readLines(Dir.CLUSTER_WORKING_DIR
                    + "/job_array/input." + jobid);

            System.out.println("starting job " + Dir.CLUSTER_WORKING_DIR
                    + "/job_array/input." + jobid + " size:" + pdbs.size());
//            logger.info("Load job input with pdb size : " + pdbs.size());

            RepeatFinder repeatFinder = null;
            //ExperRMSD finder = null;
            for (String pdb : pdbs) {
                // process and get output

                try {
                    String pdbCode = pdb.substring(0, 4);
                    String pdbChain = pdb.substring(5, 6);
                    //System.out.println(pdbCode + "\t" + pdbChain + "\n");
//                    finder = new ExperRMSD(pdbCode, pdbChain);
//                    finder.findRepeat();
//
//                    if (finder.getOutput().trim().length() > 0)
//                        buffer.append(finder.getOutput() + "\n");


                    repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                    repeatFinder.setMode(FinderMode.CONSOLE);
                    repeatFinder.findRepeat();
                    repeatFinder.setNrOfProcessors(noCore);
                    String output;

                    if (repeatFinder.isRepeat())
                        output = repeatFinder.getConsoleOutputDetails();
                    else
                        output = pdbCode + "_" + pdbChain + "\tNo-TRs";


                    buffer.append(">" + pdbCode + "_" + pdbChain + "|" + repeatFinder.getTAPOScore() + "\n" + output + "\n");

                    /**
                     * ACTIVE this code when you would like to run TAPO for large scale analysis on whole PDB
                     * To be noted that, you should switch the configuration CAL_QA within file struct-lib.cfg
                     * from on to off.
                     */
//                    output = repeatFinder.getTAPOScore();
//                    buffer.append(pdbCode + "_" + pdbChain + "\t" + repeatFinder.getTAPOScore() + "\n");


                } catch (Exception e) {
                    logger.error(e.getMessage() + " protein : " + pdb);
                }

            }

            //logger.info("job " + jobid + " finished  ");
            System.out.println("job " + jobid + " is finished");
            if (buffer.toString().length() > 0)
                DataIO.writeToFile(buffer.toString(), fileOutDir);

        } catch (Exception e) {
            logger.error(e.getMessage() + " jobid: " + jobid);
        }

    }
}
