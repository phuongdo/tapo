package apps;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.predLen.PredLenMod;
import org.cnrs.crbm.lib.repeats.RepeatRegion;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.List;
import java.util.Map;

public class ClusterJobPredLenTR {

    static Logger logger = LoggerFactory.getLogger(ClusterJobPredLenTR.class);

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

            PredLenMod predLenMod = new PredLenMod();
            // read tapo large-scale output
            ReadFasta readFasta = new ReadFasta();
            Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/resultLargeScale/PDB40_26June2015_Repeats.out");
            for (String pdb : pdbs) {
                // process and get output

                try {
                    String pdbCode = pdb.substring(0, 4);
                    String pdbChain = pdb.substring(5, 6);
                    boolean output = false;
                    if (fastaTRs.containsKey(pdb)) {
                        TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);
                        if (taPoFastaFormat.is3DRepeat()) {
                            if (taPoFastaFormat.getRepeats().size() > 0) {

                                List<RepeatRegion> regions = null;
                                try {
                                    regions = predLenMod.predictByTAPOFormat(taPoFastaFormat);
                                    for (RepeatRegion region : regions) {
                                        //Repeat repeat = taPoFastaFormat.getRepeats().get(0);
//                                        ProteinSCOP scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, region.getStart(), region.getEnd());
//                                        buffer.append(pdb + "\t" + region.getSeqStart() + "-" + region.getSeqEnd() + "\t" + scopeDomain + "\t" + region.getPredLen() + "\t" + region.getPredRU() + "\n");
//                                        output = true;
                                    }
                                } catch (Exception e) {
                                    //e.printStackTrace();
                                }

                            }
                        }
                    }
                    if (!output)
                        buffer.append(pdb + "\t" + "error\n");

                } catch (Exception e) {
//                    e.printStackTrace();
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
