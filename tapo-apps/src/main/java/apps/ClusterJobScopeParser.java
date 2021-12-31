package apps;

import org.cnrs.crbm.lib.classification.CATHParser;
import org.cnrs.crbm.lib.classification.ScopeParser;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;

public class ClusterJobScopeParser {

    static Logger logger = LoggerFactory.getLogger(ClusterJobScopeParser.class);

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
//            ScopeParser scopeParser = new ScopeParser(Dir.CLUSTER_WORKING_DIR
//                    + "/job_array/input." + jobid);
            CATHParser cathParser = new CATHParser();
            DataIO.writeToFile(cathParser.exportProperties(Dir.CLUSTER_WORKING_DIR
                    + "/job_array/input." + jobid), fileOutDir);

        } catch (Exception e) {
            e.printStackTrace();
            logger.error(e.getMessage() + " jobid: " + jobid);
        }

    }
}
