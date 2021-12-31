package org.cnrs.crbm.lib.classification;

import java.util.concurrent.ExecutionException;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.parallel.JobInput;
import org.cnrs.crbm.lib.parallel.JobOuput;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.terminal.Height;
import org.cnrs.crbm.terminal.Label;
import org.cnrs.crbm.terminal.Percentage;
import org.cnrs.crbm.terminal.ProgressBar;

/**
 * 
 * @author pdoviet
 * 
 */
public class JobDistanceCals implements Runnable {
	private int jobId;
	private JobInput jobInput;
	private Row row;

	public JobDistanceCals() {
	}

	public JobDistanceCals(int i, JobInput jobInput, Row row) {
		this.jobId = i;
		this.jobInput = jobInput;
		this.row = row;

	}

	public Object call() {

		JobOuput output = new JobOuput();
		// implement something here
		MutilAlign mutilAlign = new MutilAlign();
		int jobSize = jobInput.getProteins().size();
		int curPoint = 0;

		String threadName = Thread.currentThread().getName();

		int threadId = Integer.parseInt(threadName.split("-")[3]);
		ProgressBar progressbar = new ProgressBar(Label.create("thread:"
				+ threadId + " Job " + jobId, 25), Height.fromBottom(threadId),
				Percentage.show());
		try {
			progressbar.render().get();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		StringBuffer buffer = new StringBuffer();
		Row proteinCompare = this.row;
		try {
			Structure structure = PdbTools
					.getStructureFromLocalPdb(proteinCompare.getPdbCode());
			Atom[] atomSet1 = PdbTools.getAtomCAArray(structure
					.getChainByPDB(proteinCompare.getPdbChain()));

			for (Row protein : jobInput.getProteins()) {

				// protein name here?
				String pdbCode = protein.getPdbCode();
				String pdbChain = protein.getPdbChain();

				// do something here
				double score = 1;
				try {

					Structure structure2 = PdbTools
							.getStructureFromLocalPdb(pdbCode);
					Atom[] atomSet2 = PdbTools.getAtomCAArray(structure2
							.getChainByPDB(pdbChain));

					AFPChain afpChain = mutilAlign.pairAlign(atomSet1, atomSet2);
					// double alnLeng = (double)
					// afpChain.getOptLength();

					int sim1 = afpChain.getCoverage1();
					int sim2 = afpChain.getCoverage2();
					// score = 100 - Math.max(sim1, sim2);

					score = 1 - afpChain.getTMScore();

					buffer.append(proteinCompare.getProtein() + "\t"
							+ protein.getProtein() + "\t" + score + "\n");
					curPoint++;
					progressbar.setLabel(Label.create("thread:" + threadId
							+ " Job " + jobId + " " + pdbCode + pdbChain, 25));

					progressbar.progress((double) curPoint / jobSize).render()
							.get();

				} catch (Exception ex) {
					// ex.printStackTrace();
				}
				// output.proteins.put(pdbCode, pdbCode);

			}
		} catch (Exception ex) {

		}

		// if (buffer.toString().length() > 0)

		DataIO.writeToFile(buffer.toString(), Dir.TMP_DIR + "/output/output."
				+ jobId);

		return output;
	}

	public void run() {
		this.call();

	}

}