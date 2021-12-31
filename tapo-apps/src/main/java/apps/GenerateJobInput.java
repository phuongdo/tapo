package apps;

import com.google.common.collect.Lists;
import org.apache.commons.io.FileUtils;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;

import java.io.File;
import java.util.List;

public class GenerateJobInput {

	public static void main(String[] args) throws Exception {


        String fileInput = "data/pdbList.in";

		int partitionSize = 5;//IntMath.divide(size, nrOfJobs, RoundingMode.UP);
        //int nrOfJobs = 300;
        if (args.length > 0) {
            fileInput = args[0];
			partitionSize = Integer.parseInt(args[1]);

        }

		String dirOut = Dir.CLUSTER_WORKING_DIR + "/job_array";

		System.out.println("deleting all files and folders in  : " + dirOut);

		FileUtils.cleanDirectory(new File(dirOut));
		// ReadFasta fastaReader = new ReadFasta();
		// List<Row> rows = fastaReader.getRowAllPdbsWithChain(Dir.FASTA_LOCAL
		// + "/pdb_seqres.nrdb");

		// List<Row> rows =
		// DataIO.getListProteinFromFile("data/pdbToClassify.in");
		// List<Row> rows = DataIO.getListProteinFromFile("data/pdbList.in");

        List<String> rows = DataIO.readLines(fileInput);
        int size = rows.size();
        // int partitionSize = size / nrOfJobs;

		List<List<String>> subSets = Lists.partition(rows, partitionSize);
		System.out.println("saving " + size + " records ....");
		int jobid = 0;
		for (List<String> subs : subSets) {
			StringBuffer buffer = new StringBuffer();
			for (String row : subs) {
				buffer.append(row + "\n");
			}
			DataIO.writeToFile(buffer.toString(), dirOut + "/input." + jobid);
			jobid++;
		}

		System.out.println("done!");
	}




}
