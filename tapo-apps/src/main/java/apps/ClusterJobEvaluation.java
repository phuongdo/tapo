package apps;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.experiment.QAEvaluation;
import org.cnrs.crbm.lib.io.DataIO;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.List;

public class ClusterJobEvaluation {

    static Logger logger = LoggerFactory.getLogger(ClusterJobEvaluation.class);

    public static void main(String[] args) {

        int jobid = 0;

        if (args.length > 0)
            jobid = Integer.parseInt(args[0]);
        String fileOutDir = Dir.CLUSTER_WORKING_DIR + "/cls_output/output."
                + jobid;

        boolean isExists = new File(fileOutDir).exists();

        // System.out.println(fileOutDir);
        if (isExists) {

            System.out.println("job output " + jobid + " exists!!!");
            System.exit(1);
        }

        try {

            logger.info("Starting job : " + jobid);

            StringBuffer buffer = new StringBuffer();
            List<String> pdbs = DataIO.readLines(Dir.CLUSTER_WORKING_DIR
                    + "/job_array/input." + jobid);

            logger.info("Load job input with pdb size : " + pdbs.size());

            QAEvaluation evaluation = null;
            for (String pdb : pdbs) {
                // process and get output

                try {
                    String pdbCode = pdb.substring(0, 4);
                    String pdbChain = pdb.substring(5, 6);


                    evaluation = new QAEvaluation(pdbCode, pdbChain);
                    evaluation.findRepeat();
                    //repeatFinder.debug();
                    buffer.append(evaluation.getOutput() + "\n");


                } catch (Exception e) {
                    logger.error(e.getMessage() + " protein : " + pdb);
                }

            }

            logger.info("job " + jobid + " finished  ");

            if (buffer.toString().length() > 0)
                DataIO.writeToFile(buffer.toString(), fileOutDir);

        } catch (Exception e) {
            logger.error(e.getMessage() + " jobid: " + jobid);
        }

    }
}
