package org.cnrs.crbm.lib.classification;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import java.io.File;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

public class SeqClassify {

	public static void main(String[] args) throws Exception {

		new SeqClassify().extractSeq();

	}



	public void extractSeq() throws Exception {

		// read fasta file.
		List<Row> rows = DataIO.getListProteinFromFile("data/pdbToClassify.in");

		Set<String> pdbs = new HashSet<String>();

		for (Row row : rows) {
			pdbs.add(row.getProtein());
		}

		String dir = Dir.FASTA_LOCAL + "/pdb_seqres.txt";
		PrintWriter writer = new PrintWriter(Dir.CLUSTER_WORKING_DIR+"/pdbToCluster.fasta",
				"UTF-8");
		LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
				.readFastaProteinSequence(new File(dir));
		for (Entry<String, ProteinSequence> entry : a.entrySet()) {
			// System.out.println(entry.getValue().getOriginalHeader() + "="
			// + entry.getValue().getSequenceAsString());
			if (pdbs.contains(entry.getValue().getOriginalHeader()
					.substring(0, 6).trim()))

			{
				writer.append(">" + entry.getValue().getOriginalHeader() + "\n");
				writer.append(entry.getValue().getSequenceAsString() + "\n");

			}

		}

		writer.close();

	}
}
