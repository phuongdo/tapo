package org.cnrs.crbm.lib.parallel;

import java.util.List;
import java.util.concurrent.*;
import org.cnrs.crbm.lib.io.Row;
import com.google.common.collect.Lists;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ExecutorServiceApp {
	private int nrOfJobs = 1;
	private int nrOfProcessors = 1;
	private JobInput jobInput = new JobInput();

	public ExecutorServiceApp(int nrOfJob, int nrOfProcessors, JobInput jobInput) {
		this.nrOfJobs = nrOfJob;
		this.nrOfProcessors = nrOfProcessors;
		this.jobInput = jobInput;
	}

	// get a logger instance named "com.foo"
	//static Logger logger = Logger.getLogger(ExecutorServiceApp.class);
	static Logger logger = LoggerFactory.getLogger(ExecutorServiceApp.class);

	public JobOuput run() {

		try {
			long begTest = new java.util.Date().getTime();

			// List<Future> futuresList = new ArrayList<Future>();
			ExecutorService eservice = Executors
					.newFixedThreadPool(nrOfProcessors);

			// Divide job input to list of jobs

			int partitionSize = (int) jobInput.getProteins().size() / nrOfJobs;

			if (partitionSize == 0)
				throw new Exception(
						"number of threads must be smaller than number of threads");

			List<List<Row>> subSets = Lists.partition(jobInput.getProteins(),
					partitionSize);

			int jobid = 1;
			for (List<Row> subs : subSets) {

				JobInput partitionJobInput = new JobInput();
				partitionJobInput.setProteins(subs);
				eservice.execute(new Job(jobid, partitionJobInput));
				jobid++;
			}

			if (eservice.isTerminated()) {

				Double secs = new Double(
						(new java.util.Date().getTime() - begTest) * 0.001);

				logger.info("run time " + secs);

			}

			eservice.shutdown();

		} catch (Exception ex) {
			ex.printStackTrace();

		}
		// for (int index = 0; index < partitionJobs.size(); index++) {

		// futuresList.add(eservice.submit(new Job(index, partitionJobs
		// .get(index))));
		// }
		// get results here
		// joining output of jobs
		JobOuput joinedOutput = new JobOuput();

		// for (Future future : futuresList) {
		// try {
		// JobOuput jobOutput = (JobOuput) future.get();
		// joinedOutput.join(jobOutput);
		//
		// } catch (InterruptedException e) {
		// } catch (ExecutionException e) {
		// }
		// }

		// print the results

		// print the runtime

		return joinedOutput;

	}
}