package org.cnrs.crbm.lib.parallel;

import java.util.concurrent.ExecutionException;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.terminal.Height;
import org.cnrs.crbm.terminal.Label;
import org.cnrs.crbm.terminal.Percentage;
import org.cnrs.crbm.terminal.ProgressBar;

/**
 * 
 * @author pdoviet
 * 
 */
public class Job implements Runnable {
	private int jobId;
	private JobInput jobInput;

	public Job() {
	}

	public Job(int i, JobInput jobInput) {
		this.jobId = i;
		this.jobInput = jobInput;

	}

	public Object call() {

		JobOuput output = new JobOuput();
		// implement something here

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
		RepeatFinder repeatFinder = null;
		StringBuffer buffer = new StringBuffer();
		// buffer.append("pdbId\tfinder\tnUnits\tavgLeng\t RL\n");
		for (Row protein : jobInput.getProteins()) {

			// protein name here?
			String pdbCode = protein.getPdbCode();
			String pdbChain = protein.getPdbChain();

			try {

				// repeatFinder.debug();
				// if (repeatFinder.isRepeat() && repeatFinder.getNrMotifs() >
				// 2)
				// if (repeatFinder.isRepeat())
				// output.proteins.put(pdbCode, repeatFinder
				// .getCompareResultsBetweenEachMethos(protein));
				curPoint++;
				progressbar.setLabel(Label.create("thread:" + threadId
						+ " Job " + jobId + " " + pdbCode + pdbChain, 25));

				progressbar.progress((double) curPoint / jobSize).render()
						.get();
				repeatFinder = new RepeatFinder(pdbCode, pdbChain);
				repeatFinder.setMode(FinderMode.CONSOLE);

				repeatFinder.findRepeat();
				if (repeatFinder.isRepeat() && repeatFinder.getNrMotifs() >= 2) {

					buffer.append(repeatFinder.getConsoleOutputDetails());
				}

			} catch (Exception ex) {
				// ex.printStackTrace();
			}
			// output.proteins.put(pdbCode, pdbCode);

		}

		if (buffer.toString().length() > 0)

			DataIO.writeToFile(buffer.toString(), Dir.TMP_DIR
					+ "/output/output." + jobId);

		return output;
	}

	public void run() {
		this.call();

	}

}