package apps;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.List;

public class JobCompareStructure {
    static Logger logger = LoggerFactory.getLogger(JobCompareStructure.class);

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

            // System.out.println(Dir.FASTA_LOCAL + "/job_array/input" + jobid);
            StringBuffer buffer = new StringBuffer();
            List<String> lines = DataIO.readLines(Dir.CLUSTER_WORKING_DIR
                    + "/job_array/input." + jobid);

            logger.info("Load job input with data size : " + lines.size());

            MutilAlign mutilAlign = new MutilAlign();
            int index = 0;
            for (String line : lines) {
                // process and get output

                // System.out.print("\Rlib " + index + "/" + lines.size());
                double score = 1;
                try {
                    String pdb1 = line.split("\t")[0];
                    String pdb2 = line.split("\t")[1];

                    Atom[] atomSet1 = getAtoms(pdb1);
                    Atom[] atomSet2 = getAtoms(pdb2);
                    AFPChain afpChain = mutilAlign.pairAlignTMScoreDefault(atomSet1, atomSet2);
                    score = afpChain.getTMScore();
                    index++;
                    buffer.append(pdb1 + "\t" + pdb2 + "\t"
                            + NumberFormatUtils.format(score) + "\n");

                } catch (Exception e) {
                    logger.error(e.getMessage() + " protein : " + line);
                }

            }

            logger.info("job " + jobid + " finished  ");

            if (buffer.toString().length() > 0)
                DataIO.writeToFile(buffer.toString(), fileOutDir);

        } catch (Exception e) {
            logger.error(e.getMessage() + " jobid: " + jobid);
        }

    }


    public static Atom[] getAtoms(String pdb) {
        String pdbID = pdb.split("\\.")[0];
        String regions = pdb.split("\\.")[1];
        int start = Integer.parseInt(regions.split("_")[0]);
        int end = Integer.parseInt(regions.split("_")[1]);
        String pdbCode = pdbID.substring(0, 4);
        String pdbChain = pdbID.substring(5, 6);
        RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
        return Fragement.getFragementsofAtoms(finder.getAtoms(), start, end);
    }
}
